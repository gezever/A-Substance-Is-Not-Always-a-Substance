#!/usr/bin/env python3
"""
09f — ChemOnt equivalence assessment via Anthropic API with prompt caching.

Usage:
  python llm_chemont_equivalence.py <substances_json> <hierarchy_txt> <cache_dir> <output_json>

substances_json : JSON array of substance name strings
hierarchy_txt   : path to ChemOnt hierarchy text (used as cached system prompt)
cache_dir       : directory for per-substance JSON cache files
output_json     : path to write results array

Each result entry:
  substance_name, is_group, equivalent_node, equivalent_level,
  relationship_type, reasoning, confidence, usage
"""

import anthropic
import hashlib
import json
import os
import sys
import time


def cache_path(cache_dir: str, substance_name: str) -> str:
    key = hashlib.sha1(substance_name.encode("utf-8")).hexdigest()
    return os.path.join(cache_dir, f"{key}.json")


def load_cache(cache_dir: str, substance_name: str):
    path = cache_path(cache_dir, substance_name)
    if os.path.exists(path):
        with open(path, encoding="utf-8") as f:
            return json.load(f)
    return None


def save_cache(cache_dir: str, substance_name: str, result: dict) -> None:
    os.makedirs(cache_dir, exist_ok=True)
    with open(cache_path(cache_dir, substance_name), "w", encoding="utf-8") as f:
        json.dump(result, f, ensure_ascii=False)


PROMPT_TEMPLATE = """You are assessing whether a substance name from chemical legislation denotes a chemical class or group equivalent to a node in the ChemOnt 2.1 taxonomy.

Substance name: "{name}"

Using your chemical knowledge and the ChemOnt labels provided in the system context,
determine whether this name defines a chemical group or class (not a single specific compound), and if so, which ChemOnt node it corresponds to.

Relationship types:
- "exact"  : the name denotes the same set of chemicals as the ChemOnt class
             (e.g., "cresols" is equivalent to ChemOnt "Cresols")
- "subset" : the name denotes a specific subset of a ChemOnt class
             (e.g., "fatty acids C8-10" is a subset of ChemOnt "Fatty acids and conjugates")
- "broader": the name spans multiple ChemOnt nodes and cannot be mapped to a single one
             (e.g., "heavy metals and their compounds" covers multiple ChemOnt classes)

Return ONLY a JSON object (no markdown, no extra text):
{{
  "is_group": true or false,
  "equivalent_node": exact ChemOnt label from the list, or null,
  "equivalent_level": "kingdom" or "superclass" or "class" or "subclass" or "level 5" or "level 6", or null,
  "relationship_type": "exact" or "subset" or "broader", or null,
  "reasoning": "1-2 sentence explanation",
  "confidence": "high" or "medium" or "low"
}}

Rules:
- Use only labels that appear exactly in the ChemOnt list.
- Set is_group = false for single specific compounds, trade names, internal codes, biological extracts, and alloys.
- Set is_group = true for chemical families, element-plus-compounds groups, functional-class names, and acronym groups (PAH, PCB, PFAS, ...).
- If is_group = false, set equivalent_node, equivalent_level, and relationship_type to null.
- For "broader" relationships, set equivalent_node to the closest single ChemOnt ancestor that partially covers the group, or null if no single node is meaningful.
- Prefer a correct high-level node over a guessed specific one."""


def classify_batch(
    substances: list[str],
    hierarchy_text: str,
    cache_dir: str,
    model: str = "claude-sonnet-4-6",
    pilot: bool = False,
) -> tuple[list[dict], dict]:
    client = anthropic.Anthropic()
    results = []
    totals = {"cached": 0, "api_calls": 0,
              "input_tokens": 0, "cache_creation": 0, "cache_read": 0, "output_tokens": 0}

    for i, name in enumerate(substances):
        cached = load_cache(cache_dir, name)
        if cached is not None:
            results.append(cached)
            totals["cached"] += 1
        else:
            user_content = PROMPT_TEMPLATE.format(name=name)

            if pilot:
                print(f"\n--- PILOT substance [{i+1}]: {name!r} ---", flush=True)

            response = None
            for attempt in range(4):
                try:
                    response = client.messages.create(
                        model=model,
                        max_tokens=350,
                        system=[{
                            "type": "text",
                            "text": hierarchy_text,
                            "cache_control": {"type": "ephemeral"},
                        }],
                        messages=[{"role": "user", "content": user_content}],
                    )
                    break
                except anthropic.RateLimitError:
                    wait = 65 * (attempt + 1)
                    print(f"  Rate limit hit, waiting {wait}s (attempt {attempt+1}/3) ...", flush=True)
                    time.sleep(wait)
                    response = None
                except Exception:
                    break

            try:
                if response is None:
                    raise RuntimeError("All retry attempts failed (rate limit)")
                text = response.content[0].text.strip()
                if text.startswith("```"):
                    lines = text.splitlines()
                    text = "\n".join(lines[1:-1] if lines[-1].strip() == "```" else lines[1:])
                if pilot:
                    print(f"Response: {text}", flush=True)
                result = json.loads(text)
            except json.JSONDecodeError as e:
                result = {
                    "is_group": None,
                    "equivalent_node": None, "equivalent_level": None,
                    "relationship_type": None,
                    "reasoning": f"JSON parse error: {e}",
                    "confidence": "low",
                }
            except Exception as e:
                result = {
                    "is_group": None,
                    "equivalent_node": None, "equivalent_level": None,
                    "relationship_type": None,
                    "reasoning": f"API error: {e}",
                    "confidence": "low",
                }

            usage = {}
            if response is not None and hasattr(response, "usage"):
                u = response.usage
                usage = {
                    "input_tokens": getattr(u, "input_tokens", 0),
                    "cache_creation_input_tokens": getattr(u, "cache_creation_input_tokens", 0),
                    "cache_read_input_tokens": getattr(u, "cache_read_input_tokens", 0),
                    "output_tokens": getattr(u, "output_tokens", 0),
                }
                totals["input_tokens"]   += usage["input_tokens"]
                totals["cache_creation"] += usage["cache_creation_input_tokens"]
                totals["cache_read"]     += usage["cache_read_input_tokens"]
                totals["output_tokens"]  += usage["output_tokens"]

            result["substance_name"] = name
            result["usage"] = usage
            save_cache(cache_dir, name, result)
            results.append(result)
            totals["api_calls"] += 1

        if (i + 1) % 50 == 0 or (i + 1) == len(substances):
            print(
                f"  [{i+1}/{len(substances)}] cached={totals['cached']} "
                f"api={totals['api_calls']} "
                f"cache_read_tokens={totals['cache_read']:,}",
                flush=True,
            )

    return results, totals


def main() -> None:
    if len(sys.argv) not in (5, 6):
        print("Usage: python llm_chemont_equivalence.py "
              "<substances_json> <hierarchy_txt> <cache_dir> <output_json> [--pilot]")
        sys.exit(1)

    substances_json = sys.argv[1]
    hierarchy_txt   = sys.argv[2]
    cache_dir       = sys.argv[3]
    output_json     = sys.argv[4]
    pilot           = len(sys.argv) == 6 and sys.argv[5] == "--pilot"

    with open(substances_json, encoding="utf-8") as f:
        substances = json.load(f)

    with open(hierarchy_txt, encoding="utf-8") as f:
        hierarchy_text = f.read()

    print(f"Substances: {len(substances)}", flush=True)
    print(f"Hierarchy text: {len(hierarchy_text):,} chars (~{len(hierarchy_text)//4:,} tokens)", flush=True)
    print(f"Cache dir: {cache_dir}", flush=True)
    print(f"Model: claude-sonnet-4-6", flush=True)

    results, totals = classify_batch(substances, hierarchy_text, cache_dir, pilot=pilot)

    with open(output_json, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2, ensure_ascii=False)

    print(f"\nDone: {totals['cached']} cached, {totals['api_calls']} API calls", flush=True)
    print(f"Token usage — input: {totals['input_tokens']:,} | "
          f"cache_creation: {totals['cache_creation']:,} | "
          f"cache_read: {totals['cache_read']:,} | "
          f"output: {totals['output_tokens']:,}", flush=True)
    print(f"Results written to {output_json}", flush=True)


if __name__ == "__main__":
    main()

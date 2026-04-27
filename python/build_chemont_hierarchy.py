#!/usr/bin/env python3
"""
Build ChemOnt hierarchy text from the full SKOS TTL for use as LLM system prompt.

Usage:
  python build_chemont_hierarchy.py <skos_ttl> <output_txt>

Extracts all skos:prefLabel values and organises them by depth in the
skos:broader chain:
  depth 1  → KINGDOM
  depth 2  → SUPERCLASS
  depth 3  → CLASS
  depth 4+ → SUBCLASS

The root concept (depth 0, "Chemical entities") is excluded.
"""

import sys
import rdflib

SKOS = rdflib.Namespace("http://www.w3.org/2004/02/skos/core#")


def build_hierarchy_text(skos_ttl: str) -> str:
    g = rdflib.Graph()
    g.parse(skos_ttl, format="turtle")

    # child URI → parent URI
    broader = {str(s): str(o) for s, o in g.subject_objects(SKOS.broader)}

    # concept URI → prefLabel string
    labels = {str(s): str(o) for s, o in g.subject_objects(SKOS.prefLabel)}

    # compute depth from root (BFS-safe via memoisation)
    memo = {}

    def depth(uri):
        if uri in memo:
            return memo[uri]
        stack, node = [], uri
        while node not in memo:
            stack.append(node)
            if node not in broader:
                memo[node] = 0
                break
            node = broader[node]
        for n in reversed(stack):
            memo[n] = 1 + memo.get(broader.get(n, n), 0) if n in broader else 0
        return memo[uri]

    max_depth = max((depth(uri) for uri in labels if depth(uri) > 0), default=11)
    by_level = {d: [] for d in range(1, max_depth + 1)}
    for uri, label in labels.items():
        d = depth(uri)
        if d == 0:
            continue  # skip root "Chemical entities"
        by_level[d].append(label)

    for d in by_level:
        by_level[d] = sorted(set(by_level[d]))

    def level_name(d):
        names = {1: "KINGDOM", 2: "SUPERCLASS", 3: "CLASS", 4: "SUBCLASS"}
        return names.get(d, f"LEVEL {d}")

    lines = [
        "=== ChemOnt 2.1 Chemical Taxonomy (Full) ===",
        "",
        "Use ONLY the exact labels listed below.",
        "Taxonomy levels (broadest to most specific):",
        "  kingdom → superclass → class → subclass → level 5 → level 6 → ...",
        "Use null when uncertain rather than guessing a wrong label.",
    ]

    for d in sorted(by_level):
        if not by_level[d]:
            continue
        lines += ["", f"--- {level_name(d)} ---"]
        lines += [f"  {lbl}" for lbl in by_level[d]]

    text = "\n".join(lines)

    counts = {level_name(b): len(by_level[b]) for b in by_level}
    total = sum(counts.values())
    print("Labels per level:")
    for name, n in counts.items():
        print(f"  {name}: {n}")
    print(f"  Total: {total}")
    print(f"Text length: {len(text):,} chars (~{len(text)//4:,} tokens estimated)")

    return text


def main():
    if len(sys.argv) != 3:
        print("Usage: python build_chemont_hierarchy.py <skos_ttl> <output_txt>")
        sys.exit(1)

    skos_ttl, output_txt = sys.argv[1], sys.argv[2]
    text = build_hierarchy_text(skos_ttl)

    with open(output_txt, "w", encoding="utf-8") as f:
        f.write(text)

    print(f"Written to {output_txt}")


if __name__ == "__main__":
    main()

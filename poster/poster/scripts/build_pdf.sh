#!/usr/bin/env bash
set -euo pipefail


#
# cp ../../output/figures/Analysis_4bd_v_Entity_type_and_linkability_vertical.pdf ../images/.

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TEX_FILE="${1:-poster.tex}"

cd "$ROOT_DIR"

if [[ ! -f "$TEX_FILE" ]]; then
  echo "Error: TeX file '$TEX_FILE' not found." >&2
  exit 1
fi

BASE_NAME="${TEX_FILE%.tex}"

# poster.tex requires LuaTeX (fontspec/Flanders Art fonts); minted requires
# shell-escape; latexmk handles bibtex/passes automatically.
if command -v latexmk >/dev/null 2>&1; then
  latexmk -pdflua -shell-escape -interaction=nonstopmode -file-line-error "$TEX_FILE"
  echo "Build finished: ${BASE_NAME}.pdf"
  cp "${BASE_NAME}.pdf" poster_24.pdf
  exit 0
fi

if ! command -v lualatex >/dev/null 2>&1; then
  echo "Error: neither 'latexmk' nor 'lualatex' was found in PATH." >&2
  exit 1
fi

if ! command -v bibtex >/dev/null 2>&1; then
  echo "Error: 'bibtex' was not found in PATH." >&2
  exit 1
fi

lualatex -shell-escape -interaction=nonstopmode -file-line-error "$TEX_FILE"

if rg -n "\\\\bibliography\s*\{" "$TEX_FILE" >/dev/null 2>&1; then
  bibtex "$BASE_NAME"
fi

lualatex -shell-escape -interaction=nonstopmode -file-line-error "$TEX_FILE"
lualatex -shell-escape -interaction=nonstopmode -file-line-error "$TEX_FILE"

echo "Build finished: ${BASE_NAME}.pdf"



#!/usr/bin/env bash
set -euo pipefail


cp ../../output/figures/Analysis_4bd_v_Entity_type_and_linkability_vertical.pdf ../images/.
cp ../../output/figures/Analysis_10c_workload_non_structured.pdf ../images/.
cp ../../output/figures/Analysis_14e_Bubble_priority.pdf ../images/.
cp ../../output/figures/Analysis_89_UMAP_and_ChemOnt_quality_ellipse.pdf ../images/.
cp ../../output/figures/Analysis_10c_workload_curve_per_source.pdf ../images/.

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TEX_FILE="${1:-main.tex}"

cd "$ROOT_DIR"

if [[ ! -f "$TEX_FILE" ]]; then
  echo "Error: TeX file '$TEX_FILE' not found." >&2
  exit 1
fi

BASE_NAME="${TEX_FILE%.tex}"

# minted requires shell-escape; latexmk handles bibtex/passes automatically.
if command -v latexmk >/dev/null 2>&1; then
  latexmk -pdf -shell-escape -interaction=nonstopmode -file-line-error "$TEX_FILE"
  exit 0
fi

if ! command -v pdflatex >/dev/null 2>&1; then
  echo "Error: neither 'latexmk' nor 'pdflatex' was found in PATH." >&2
  exit 1
fi

if ! command -v bibtex >/dev/null 2>&1; then
  echo "Error: 'bibtex' was not found in PATH." >&2
  exit 1
fi

pdflatex -shell-escape -interaction=nonstopmode -file-line-error "$TEX_FILE"

if rg -n "\\\\bibliography\s*\{" "$TEX_FILE" >/dev/null 2>&1; then
  bibtex "$BASE_NAME"
fi

pdflatex -shell-escape -interaction=nonstopmode -file-line-error "$TEX_FILE"
pdflatex -shell-escape -interaction=nonstopmode -file-line-error "$TEX_FILE"

echo "Build finished: ${BASE_NAME}.pdf"

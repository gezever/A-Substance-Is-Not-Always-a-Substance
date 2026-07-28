#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TEX_FILE="${1:-main.tex}"
FULL_CLEAN="${2:-}"

cd "$ROOT_DIR"

if [[ ! -f "$TEX_FILE" ]]; then
  echo "Error: TeX file '$TEX_FILE' not found." >&2
  exit 1
fi

if ! command -v latexmk >/dev/null 2>&1; then
  echo "Error: 'latexmk' was not found in PATH (required for cleanup script)." >&2
  exit 1
fi

if [[ "$FULL_CLEAN" == "--all" ]]; then
  latexmk -C "$TEX_FILE"
  echo "Removed auxiliary files and generated PDF for $TEX_FILE"
else
  latexmk -c "$TEX_FILE"
  echo "Removed auxiliary files for $TEX_FILE"
fi

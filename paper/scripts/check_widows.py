#!/usr/bin/env python3
"""
Tool to check for "widow lines" in the PDF output (lines with very few characters)
(c) 2022 Maxim Van de Wynckel - Vrije Universiteit Brussel
"""
import re
import sys
from collections import defaultdict

try:
    from pypdf import PdfReader
except ImportError:
    print("pypdf not installed. Run: pip3 install pypdf")
    sys.exit(1)

## ---------------------------- ##
## CONFIGURE YOUR PAPER HERE
PDF = "main.pdf"
MAX_BODY_PAGES = 15        # ignore bibliography-heavy tail
MAX_CHARS = 40             # flag lines shorter than this
MIN_GAP_PT = 8             # minimum gap to next line (filters mid-paragraph wraps)
MIN_FONTSIZE = 7           # ignore footnote/caption small text
MAX_FONTSIZE = 11.5        # ignore headings (larger font)
## ---------------------------- ##

IGNORE_PATTERNS = [
    r"^\d+\.",                              # section headings "1. Introduction"
    r"^[\d\s\.,\-\+×]+$",                   # axis numbers / tick labels
    r"^\(A\)|\(B\)|\(C\)",                  # figure panel labels
    r"^(References|Data Av|Suppl|LLM)",     # back-matter headings
    r"^(Model|Hit@|MRR)",                   # table headers
]

reader = PdfReader(PDF)
issues = []

for page_num, page in enumerate(reader.pages, 1):
    if page_num > MAX_BODY_PAGES:
        break

    tokens = []

    def visitor(text, cm, tm, fontdict, fontsize):
        if text.strip():
            tokens.append((tm[5], tm[4], text, fontsize or 0))

    page.extract_text(visitor_text=visitor)
    if not tokens:
        continue

    lines = defaultdict(list)
    for y, x, text, fs in tokens:
        lines[round(y / 2) * 2].append((x, text, fs))

    sorted_lines = sorted(lines.items(), key=lambda t: -t[0])

    for i, (y, chunks) in enumerate(sorted_lines[:-1]):
        chunks_s = sorted(chunks, key=lambda c: c[0])
        line_text = " ".join(c[1] for c in chunks_s).strip()
        avg_fs = sum(c[2] for c in chunks_s) / len(chunks_s)

        if not (MIN_FONTSIZE <= avg_fs <= MAX_FONTSIZE):
            continue
        if any(re.search(p, line_text) for p in IGNORE_PATTERNS):
            continue
        if len(line_text) > MAX_CHARS:
            continue

        gap_below = y - sorted_lines[i + 1][0]
        if gap_below < MIN_GAP_PT:
            continue

        severity = "CRITICAL" if len(line_text) <= 15 else "WARNING"
        issues.append((severity, page_num, len(line_text), gap_below, line_text))

print(f"=== Widow line report for {PDF} ===\n")
if not issues:
    print("No widow lines detected.")
else:
    for sev, pn, ln, gap, text in sorted(issues, key=lambda x: (x[1], x[2])):
        tag = f"[{sev}]"
        print(f"  p{pn:2d} {tag:10s} [{ln:2d} chars, gap={gap:.0f}pt]: '{text}'")
    print(f"\nTotal: {len(issues)}  "
          f"({sum(1 for i in issues if i[0]=='CRITICAL')} critical, "
          f"{sum(1 for i in issues if i[0]=='WARNING')} warnings)")

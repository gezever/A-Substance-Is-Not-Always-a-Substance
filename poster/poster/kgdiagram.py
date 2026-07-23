#!/usr/bin/env python3
"""Build the central knowledge-graph diagram as TikZ, with layout validation.
Coordinates are in mm on a fixed canvas; everything is checked for overlap
and containment before a single line of TikZ is written."""

import itertools
import math

CANVAS_W, CANVAS_H = 780.0, 312.0

# style -> (fill, draw, dashed, bold title)
STYLES = {
    "source": ("volight", "vodark", False),
    "hub":    ("voyellow", "vodark", False),
    "layer":  ("white",   "vodark", False),
    "soft":   ("volight", "vodark", False),
    "sep":    ("white",   "vodark", True),
    "legend": ("none",    "vodark", False),
}

# name: (x, y, w, h, style, title, detail lines)
BOXES = {
    # ---- tier A : external sources -------------------------------------
    "A1": (2.5, 255, 147, 55, "source", "18 regulatory lists",
           ["ECHA · Pesticides", "OSPAR · VMM"]),
    "A2": (159.5, 255, 147, 55, "source", "PubChem",
           ["CAS \\textrightarrow{} InChIKey"]),
    "A3": (316.5, 255, 147, 55, "source", "ClassyFire",
           ["InChIKey \\textrightarrow{}", "ChemOnt direct parent"]),
    "A4": (473.5, 255, 147, 55, "source", "ChEBI ontology",
           ["EBI · OWL"]),
    "A5": (630.5, 255, 147, 55, "source", "Embeddings + LLM",
           ["MiniLM · $k$-means", "Claude Sonnet 4.6"]),

    # ---- tier B : the hub ----------------------------------------------
    "B0": (2.5, 175, 165, 55, "layer", "skos:ConceptScheme",
           ["chemical\\_substance"]),
    "HUB": (192, 170, 396, 65, "hub", "Regulatory substances",
            ["${\\sim}$17,500 skos:Concept + dbo:ChemicalSubstance",
             "CAS · EC number · InChIKey"]),
    "SEP": (612.5, 170, 165, 65, "sep", "Scope mapping",
            ["304 SKOS triples", "kept separate,", "pending validation"]),

    # ---- tier C1 : semantic layers -------------------------------------
    "C1": (2.5, 86, 180, 62, "soft", "13 skos:Collection",
           ["7 linkability tiers", "6 embedding clusters"]),
    "C2": (200, 86, 180, 62, "layer", "ChemOnt 2.1",
           ["${\\sim}$4,800 skos:Concept", "Kingdom \\textrightarrow{} SubClass"]),
    "C3": (397.5, 86, 180, 62, "layer", "ChEBI substances",
           ["${\\sim}$4,000 owl:Class", "SMILES · InChI · mass"]),
    "C4": (595, 86, 182, 62, "layer", "ChEBI biological roles",
           ["${\\sim}$1,400 owl:Class", "carcinogen · allergen"]),

    # ---- tier C2 -------------------------------------------------------
    "D2": (200, 8, 180, 58, "soft", "xkos:ClassificationLevel",
           ["Kingdom · SuperClass", "Class · SubClass"]),
    "D3": (397.5, 8, 379.5, 58, "soft", "oa:Annotation",
           ["${\\sim}$46,600 cross-domain links",
            "regulatory provenance + biological roles"]),
}

# (from, to, label, side)  side: 'l'/'r'/'c' nudges the label off the line.
# Routing (straight vs. bent around an intervening tier) is computed, not
# hand-placed -- see route() below.
EDGES = [
    ("A1", "HUB", "import + normalisation", "l"),
    ("A2", "HUB", "resolve identity", "l"),
    ("A3", "C2",  "skos:broader", "r"),
    ("A4", "C3",  "riot OWL\\textrightarrow{}TTL", "r"),
    ("A5", "C1",  "clusters", "r"),
    ("A5", "SEP", "", "c"),
    ("HUB", "B0", "skos:inScheme", "c"),
    ("HUB", "C2", "skos:broader", "l"),
    ("HUB", "C3", "skos:exactMatch", "r"),
    ("C1", "HUB", "skos:member", "l"),
    ("C2", "D2",  "xkos:inLevel", "c"),
    ("C3", "C4",  "rdfs:subClassOf", "c"),
    ("D3", "C4",  "oa:hasBody", "c"),
    ("HUB", "SEP", "", "c"),
]

# ------------------------------------------------------------ tier layout
# Rows are laid out left-to-right independently; D2/D3 are not part of the
# search because their geometry is derived from C2/C3/C4 after the fact
# (D2 sits under whichever box ends up in C2's slot; D3 spans C3+C4).
ROW_A = ["A1", "A2", "A3", "A4", "A5"]
ROW_B = ["B0", "HUB", "SEP"]
ROW_C = ["C1", "C2", "C3", "C4"]
ROWS = [ROW_A, ROW_B, ROW_C, ["D2", "D3"]]
ROW_OF = {name: i for i, row in enumerate(ROWS) for name in row}

def row_gap(row):
    """The (uniform, by construction) gap between consecutive boxes in a row."""
    ordered = sorted(row, key=lambda n: BOXES[n][0])
    n1, n2 = ordered[0], ordered[1]
    return BOXES[n2][0] - (BOXES[n1][0] + BOXES[n1][2])

def layout_row(order, start_x, gap):
    """x-position for each box name in `order`, packed left-to-right from
    start_x using each box's own width; total row width is order-independent
    since it only depends on the (fixed) set of widths and the gap."""
    xs, cur = {}, start_x
    for n in order:
        xs[n] = cur
        cur += BOXES[n][2] + gap
    return xs

def total_edge_length(xs):
    def centre(n):
        x = xs.get(n, BOXES[n][0])
        return (x + BOXES[n][2] / 2, BOXES[n][1] + BOXES[n][3] / 2)
    total = 0.0
    for a, b, _lab, _side in EDGES:
        (ax, ay), (bx, by) = centre(a), centre(b)
        total += math.hypot(ax - bx, ay - by)
    return total

def optimize_order():
    """Brute-force the left-right order within each of ROW_A/B/C (independent
    permutations, small enough to enumerate fully) to minimise total edge
    length. C3/C4 are kept adjacent so D3 can still span exactly those two
    columns. Returns {name: new_x} for every box in ROW_A + ROW_B + ROW_C."""
    start_a, gap_a = BOXES["A1"][0], row_gap(ROW_A)
    start_b, gap_b = BOXES["B0"][0], row_gap(ROW_B)
    start_c, gap_c = BOXES["C1"][0], row_gap(ROW_C)

    best_cost, best_xs = None, None
    for perm_a in itertools.permutations(ROW_A):
        xs_a = layout_row(perm_a, start_a, gap_a)
        for perm_b in itertools.permutations(ROW_B):
            xs_b = layout_row(perm_b, start_b, gap_b)
            for perm_c in itertools.permutations(ROW_C):
                if abs(perm_c.index("C3") - perm_c.index("C4")) != 1:
                    continue
                xs_c = layout_row(perm_c, start_c, gap_c)
                xs = {**xs_a, **xs_b, **xs_c}
                cost = total_edge_length(xs)
                if best_cost is None or cost < best_cost:
                    best_cost, best_xs = cost, xs
    return best_xs

def apply_order(xs):
    for name, x in xs.items():
        old = list(BOXES[name])
        old[0] = x
        BOXES[name] = tuple(old)
    # D2 tracks C2's slot exactly (same x and width).
    c2 = BOXES["C2"]
    d2 = list(BOXES["D2"]); d2[0], d2[2] = c2[0], c2[2]
    BOXES["D2"] = tuple(d2)
    # D3 spans the bounding box of C3 and C4 (kept adjacent by construction).
    c3, c4 = BOXES["C3"], BOXES["C4"]
    left = min(c3[0], c4[0])
    right = max(c3[0] + c3[2], c4[0] + c4[2])
    d3 = list(BOXES["D3"]); d3[0], d3[2] = left, right - left
    BOXES["D3"] = tuple(d3)

apply_order(optimize_order())

# ------------------------------------------------------------------ checks
def rect(b):
    x, y, w, h = b[0], b[1], b[2], b[3]
    return (x, y, x + w, y + h)

def overlaps(r1, r2, pad=0.0):
    return not (r1[2] + pad <= r2[0] or r2[2] + pad <= r1[0]
                or r1[3] + pad <= r2[1] or r2[3] + pad <= r1[1])

def validate():
    errs = []
    names = list(BOXES)
    for i, n in enumerate(names):
        r = rect(BOXES[n])
        if r[0] < 0 or r[1] < 0 or r[2] > CANVAS_W or r[3] > CANVAS_H:
            errs.append(f"{n} outside canvas: {r}")
        for m in names[i+1:]:
            if overlaps(r, rect(BOXES[m]), pad=6.0):
                errs.append(f"{n} and {m} overlap (or <6mm apart)")
    return errs

# ------------------------------------------------------------------ emit
def anchors(a, b):
    """Pick sensible exit/entry points on the two box borders."""
    ax, ay, aw, ah = BOXES[a][:4]
    bx, by, bw, bh = BOXES[b][:4]
    acx, acy = ax + aw/2, ay + ah/2
    bcx, bcy = bx + bw/2, by + bh/2
    dx, dy = bcx - acx, bcy - acy
    if abs(dy) >= abs(dx):                     # mostly vertical
        if dy < 0:
            p1 = (acx, ay); p2 = (bcx, by + bh)
        else:
            p1 = (acx, ay + ah); p2 = (bcx, by)
    else:                                      # mostly horizontal
        if dx > 0:
            p1 = (ax + aw, acy); p2 = (bx, bcy)
        else:
            p1 = (ax, acy); p2 = (bx + bw, bcy)
    return p1, p2

def nearest_gap(row, near_x):
    """x-coordinate of the midpoint of the row's gap closest to near_x."""
    ordered = sorted(row, key=lambda n: BOXES[n][0])
    gaps = []
    for n1, n2 in zip(ordered, ordered[1:]):
        left = BOXES[n1][0] + BOXES[n1][2]
        right = BOXES[n2][0]
        gaps.append((left + right) / 2)
    return min(gaps, key=lambda g: abs(g - near_x))

def route(a, b, stagger=0.0):
    """Path (list of points) for edge a->b. Adjacent/same-tier edges get a
    straight line. An edge that skips a tier (currently: tier A -> tier C)
    is bent through the nearest free gap in the tier it would otherwise cut
    across, so it never crosses a box that sits between its endpoints."""
    p1, p2 = anchors(a, b)
    if abs(ROW_OF[a] - ROW_OF[b]) < 2:
        return [p1, p2]
    mid_row = ROWS[(ROW_OF[a] + ROW_OF[b]) // 2]
    gap_x = nearest_gap(mid_row, (p1[0] + p2[0]) / 2) + stagger
    band_hi = (min(BOXES[n][1] for n in ROWS[ROW_OF[a]])
               + max(BOXES[n][1] + BOXES[n][3] for n in mid_row)) / 2 + stagger
    top_b = max(BOXES[n][1] + BOXES[n][3] for n in ROWS[ROW_OF[b]])
    return [p1, (p1[0], band_hi), (gap_x, band_hi), (gap_x, top_b), p2]

def render():
    errs = validate()
    if errs:
        raise SystemExit("LAYOUT ERRORS:\n  " + "\n  ".join(errs))

    L = []
    L.append("%% central knowledge-graph diagram (generated by kgdiagram.py)")
    L.append("\\begin{tikzpicture}[x=1mm,y=1mm]")
    L.append(f"\\path[use as bounding box] (0,0) rectangle ({CANVAS_W:.1f},{CANVAS_H:.1f});")

    # boxes first (bottom layer), so lines and labels always paint on top of them
    for name, (x, y, w, h, style, title, details) in BOXES.items():
        fill, draw, dashed = STYLES[style]
        opts = [f"anchor=south west", f"minimum width={w}mm", f"minimum height={h}mm",
                f"text width={w-12}mm", "align=center", "inner sep=0pt",
                f"draw={draw}", "line width=0.8mm", "rounded corners=0pt"]
        if fill != "none":
            opts.append(f"fill={fill}")
        if dashed:
            opts.append("dash pattern=on 5mm off 3.5mm")
        body = ("\\fontsize{25}{28}\\selectfont\\bfseries " + title)
        if details:
            body += ("\\\\[1.2mm]\\fontsize{19}{23}\\selectfont\\mdseries "
                     + "\\\\ ".join(details))
        L.append(f"\\node[{', '.join(opts)}] at ({x},{y}) {{{body}}};")

    # edges on top of boxes, so a line is never hidden by a box it happens to pass near.
    # Edges that skip a tier are staggered (fanned out around the direct path)
    # so their bent routes don't run on top of one another.
    skip_edges = [i for i, (a, b, *_r) in enumerate(EDGES) if abs(ROW_OF[a] - ROW_OF[b]) >= 2]
    stagger_of = {i: (k - (len(skip_edges) - 1) / 2) * 6.0 for k, i in enumerate(skip_edges)}

    edge_labels = []
    for i, (a, b, lab, side) in enumerate(EDGES):
        pts = route(a, b, stagger=stagger_of.get(i, 0.0))
        dash = "dash pattern=on 4mm off 3mm, " if (a == "SEP" or b == "SEP") else ""
        path = " -- ".join(f"({x:.1f},{y:.1f})" for x, y in pts)
        L.append(f"\\draw[{dash}-{{Stealth[length=5mm,width=3.5mm]}}, line width=0.9mm, vodark] "
                 f"{path};")
        if lab:
            # place the label on the longest segment of the (possibly bent) path
            segs = list(zip(pts, pts[1:]))
            (sx1, sy1), (sx2, sy2) = max(segs, key=lambda s: (s[0][0]-s[1][0])**2 + (s[0][1]-s[1][1])**2)
            mx, my = (sx1+sx2)/2, (sy1+sy2)/2
            off = {"l": (-3, 0), "r": (3, 0), "c": (0, 0)}[side]
            anch = {"l": "east", "r": "west", "c": "center"}[side]
            edge_labels.append((mx + off[0], my + off[1], anch, lab))

    # labels last (top layer), so a label is always legible over both boxes and lines
    for mx, my, anch, lab in edge_labels:
        L.append(f"\\node[anchor={anch}, inner sep=1.2mm, fill=white, draw=vodark, line width=0.3mm, "
                 f"font=\\fontsize{{17}}{{20}}\\selectfont] "
                 f"at ({mx:.1f},{my:.1f}) {{\\ttfamily {lab}}};")

    L.append("\\end{tikzpicture}%")
    return "\n".join(L) + "\n"

def emit(path="kgdiagram.tex"):
    text = render()
    open(path, "w").write(text)
    return len(BOXES), len(EDGES)

if __name__ == "__main__":
    import sys
    if "--print" in sys.argv:
        sys.stdout.write(render())
        raise SystemExit(0)
    nb, ne = emit()
    print(f"layout valid: {nb} boxes, {ne} edges, canvas {CANVAS_W}x{CANVAS_H} mm")
    used = max(rect(b)[3] for b in BOXES.values())
    print(f"top-most content: {used} mm ; right-most: {max(rect(b)[2] for b in BOXES.values()):.1f} mm")

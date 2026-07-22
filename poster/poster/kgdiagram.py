#!/usr/bin/env python3
"""Build the central knowledge-graph diagram as TikZ, with layout validation.
Coordinates are in mm on a fixed canvas; everything is checked for overlap
and containment before a single line of TikZ is written."""

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

# (from, to, label, side)  side: 'l'/'r'/'c' nudges the label off the line
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

def emit(path="kgdiagram.tex"):
    errs = validate()
    if errs:
        raise SystemExit("LAYOUT ERRORS:\n  " + "\n  ".join(errs))

    L = []
    L.append("%% central knowledge-graph diagram (generated by kgdiagram.py)")
    L.append("\\begin{tikzpicture}[x=1mm,y=1mm]")
    L.append(f"\\path[use as bounding box] (0,0) rectangle ({CANVAS_W:.1f},{CANVAS_H:.1f});")

    # edges first so boxes paint over the ends
    for a, b, lab, side in EDGES:
        p1, p2 = anchors(a, b)
        dash = "dash pattern=on 4mm off 3mm, " if (a == "SEP" or b == "SEP") else ""
        L.append(f"\\draw[{dash}-{{Stealth[length=5mm,width=3.5mm]}}, line width=0.9mm, vodark] "
                 f"({p1[0]:.1f},{p1[1]:.1f}) -- ({p2[0]:.1f},{p2[1]:.1f});")
        if lab:
            mx, my = (p1[0]+p2[0])/2, (p1[1]+p2[1])/2
            off = {"l": (-3, 0), "r": (3, 0), "c": (0, 0)}[side]
            anch = {"l": "east", "r": "west", "c": "center"}[side]
            L.append(f"\\node[anchor={anch}, inner sep=1.2mm, fill=white, "
                     f"font=\\fontsize{{17}}{{20}}\\selectfont] "
                     f"at ({mx+off[0]:.1f},{my+off[1]:.1f}) {{\\ttfamily {lab}}};")

    # boxes
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

    L.append("\\end{tikzpicture}%")
    open(path, "w").write("\n".join(L) + "\n")
    return len(BOXES), len(EDGES)

if __name__ == "__main__":
    nb, ne = emit()
    print(f"layout valid: {nb} boxes, {ne} edges, canvas {CANVAS_W}x{CANVAS_H} mm")
    used = max(rect(b)[3] for b in BOXES.values())
    print(f"top-most content: {used} mm ; right-most: {max(rect(b)[2] for b in BOXES.values()):.1f} mm")

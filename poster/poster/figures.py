#!/usr/bin/env python3
"""Redraw figures 1 and 2 as TikZ (vector), replacing the raster crops.
Layout is computed and validated in mm before any TikZ is emitted."""

# --------------------------------------------------------------- shared
FILLS = {
    "chem":    "volight",
    "reg":     "voyellow!35",
    "example": "volight",
    "plain":   "white",
    "sum":     "voyellow!35",
    "box":     "white",
}

def rects_overlap(a, b, pad=0.0):
    return not (a[0]+a[2]+pad <= b[0] or b[0]+b[2]+pad <= a[0]
                or a[1]+a[3]+pad <= b[1] or b[1]+b[3]+pad <= a[1])


class Fig:
    def __init__(self, w, h):
        self.W, self.H = w, h
        self.nodes = {}      # name -> (x, y, w, h, style, lines, fontsize)
        self.groups = {}     # name -> (x, y, w, h, label)
        self.edges = []      # (from, to, kind)

    def node(self, name, x, y, w, h, style, lines, fs=15):
        self.nodes[name] = (x, y, w, h, style, lines, fs)

    def group(self, name, x, y, w, h, label):
        self.groups[name] = (x, y, w, h, label)

    def edge(self, a, b, kind="elbow"):
        self.edges.append((a, b, kind))

    def validate(self):
        errs = []
        for n, v in self.nodes.items():
            x, y, w, h = v[:4]
            if x < 0 or y < 0 or x+w > self.W or y+h > self.H:
                errs.append(f"node {n} outside canvas")
        names = list(self.nodes)
        for i, a in enumerate(names):
            for b in names[i+1:]:
                if rects_overlap(self.nodes[a][:4], self.nodes[b][:4], pad=2.0):
                    errs.append(f"nodes {a} and {b} overlap")
        for g, v in self.groups.items():
            x, y, w, h = v[:4]
            if x < 0 or y < 0 or x+w > self.W or y+h > self.H:
                errs.append(f"group {g} outside canvas")
        return errs

    # ---- geometry helpers
    def cx(self, n): return self.nodes[n][0] + self.nodes[n][2]/2
    def cy(self, n): return self.nodes[n][1] + self.nodes[n][3]/2
    def right(self, n): return self.nodes[n][0] + self.nodes[n][2]
    def left(self, n): return self.nodes[n][0]
    def top(self, n): return self.nodes[n][1] + self.nodes[n][3]
    def bot(self, n): return self.nodes[n][1]

    def emit(self, path):
        errs = self.validate()
        if errs:
            raise SystemExit("LAYOUT ERRORS:\n  " + "\n  ".join(errs))
        L = ["\\begin{tikzpicture}[x=1mm,y=1mm]",
             f"\\path[use as bounding box] (0,0) rectangle ({self.W},{self.H});"]

        # group containers first (behind)
        for g, (x, y, w, h, label) in self.groups.items():
            L.append(f"\\draw[fill=white, draw=vodark, line width=0.5mm, rounded corners=1mm] "
                     f"({x},{y}) rectangle ({x+w},{y+h});")
            L.append(f"\\node[anchor=north west, inner sep=0pt, font=\\fontsize{{16}}{{19}}"
                     f"\\selectfont\\bfseries, text=vodark] at ({x+4},{y+h-3}) {{{label}}};")

        # edges
        for a, b, kind in self.edges:
            if kind == "elbow":                    # left -> right tree branch
                x1, y1 = self.right(a), self.cy(a)
                x2, y2 = self.left(b), self.cy(b)
                xm = (x1 + x2) / 2
                pts = f"({x1},{y1}) -- ({xm},{y1}) -- ({xm},{y2}) -- ({x2},{y2})"
            elif kind == "down":                   # top -> bottom, vertical elbow
                x1, y1 = self.cx(a), self.bot(a)
                x2, y2 = self.cx(b), self.top(b)
                ym = (y1 + y2) / 2
                pts = f"({x1},{y1}) -- ({x1},{ym}) -- ({x2},{ym}) -- ({x2},{y2})"
            else:                                  # straight
                x1, y1 = self.right(a), self.cy(a)
                x2, y2 = self.left(b), self.cy(b)
                pts = f"({x1},{y1}) -- ({x2},{y2})"
            L.append(f"\\draw[-{{Stealth[length=3.2mm,width=2.4mm]}}, line width=0.5mm, "
                     f"vodark, rounded corners=1.5mm] {pts};")

        # nodes on top
        for n, (x, y, w, h, style, lines, fs) in self.nodes.items():
            body = ("\\\\ ".join(lines))
            L.append(f"\\node[anchor=south west, minimum width={w}mm, minimum height={h}mm, "
                     f"text width={w-6}mm, align=center, inner sep=0pt, "
                     f"fill={FILLS[style]}, draw=vodark, line width=0.5mm, rounded corners=1mm, "
                     f"font=\\fontsize{{{fs}}}{{{int(fs*1.2)}}}\\selectfont, text=vodark] "
                     f"at ({x},{y}) {{{body}}};")
        L.append("\\end{tikzpicture}%")
        open(path, "w").write("\n".join(L) + "\n")
        return len(self.nodes), len(self.edges)


# ==================================================================== FIGURE 1
def figure1():
    f = Fig(246, 110)
    # ---- regulatory container (left), 3 columns x 4 rows
    f.group("reg", 0, 0, 168, 102, "Regulatory datasets")
    rows_y = [6, 27, 48, 69]                       # bottom-up: H, G, F, E
    tw, ew = 46, 60
    tx, ex = 46, 100
    f.node("D", 5, 38, 34, 18, "reg", ["Substance"], 16)
    specs = [("H", "Analytical parameter", "H1", ["Example:", "Sum PFOS + PFBA"]),
             ("G", "Substance group",      "G1", ["Example:", "Perfluorinated compounds"]),
             ("F", "Mixture",              "F1", ["Example:", "commercial PFAS mixture"]),
             ("E", "Individual molecule",  "E1", ["Example: PFOS"])]
    for (tn, tlabel, en, elines), y in zip(specs, rows_y):
        f.node(tn, tx, y, tw, 18, "reg", [tlabel], 14)
        f.node(en, ex, y, ew, 18, "example", elines, 13)
        f.edge("D", tn, "elbow")
        f.edge(tn, en, "elbow")

    # ---- chemistry container (right), 3 stacked boxes
    f.group("chem", 176, 0, 70, 102, "Chemistry")
    f.node("C", 181, 6,  60, 22, "chem", ["Structure identifiers", "InChI / InChIKey"], 13)
    f.node("B", 181, 34, 60, 22, "chem", ["Structure representation", "atoms + bonds"], 13)
    f.node("A", 181, 62, 60, 18, "chem", ["Molecule"], 15)
    f.edge("A", "B", "down")
    f.edge("B", "C", "down")
    return f


# ==================================================================== FIGURE 2
def figure2():
    f = Fig(246, 116)
    f.node("U", 63, 94, 120, 18, "plain", ["Universe of chemical molecules"], 15)

    f.group("pfas", 0, 38, 246, 46, "Regulatory group: PFAS")
    mw, gap = 56.5, 5
    xs = [5 + i*(mw+gap) for i in range(4)]
    labels = [["PFOS molecule"], ["PFBA molecule"], ["PFBS molecule"],
              ["Other fluorinated", "molecules"]]
    names = ["PFOS", "PFBA", "PFBS", "OTHER"]
    for n, x, lab in zip(names, xs, labels):
        f.node(n, x, 43, mw, 24, "chem", lab, 14)
        f.edge("U", n, "down")

    f.group("ana", 0, 0, 140, 32, "Analytical parameter")
    f.node("SUM", 5, 4, 100, 18, "sum", ["Sum PFOS + PFBA"], 15)
    f.edge("PFOS", "SUM", "down")
    f.edge("PFBA", "SUM", "down")
    return f


if __name__ == "__main__":
    for name, fn in [("fig1.tex", figure1), ("fig2.tex", figure2)]:
        fig = fn()
        n, e = fig.emit(name)
        print(f"{name}: layout valid, {n} nodes, {e} edges, "
              f"canvas {fig.W}x{fig.H} mm")

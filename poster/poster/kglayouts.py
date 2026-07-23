#!/usr/bin/env python3
"""Three alternative layouts for the central knowledge-graph diagram.

Every layout is checked before emission:
  * no node outside the canvas
  * no two nodes overlapping
  * no edge segment passing through a node
  * edge/edge crossings counted and reported
Edges are routed orthogonally through reserved channels.
"""

# ------------------------------------------------------------------ content
# key: (title, [detail lines])
TEXT = {
    "SRC_REG":   ("18 regulatory lists", ["ECHA · Pesticides", "OSPAR · VMM"]),
    "SRC_PUB":   ("PubChem",             ["CAS \\textrightarrow{} InChIKey"]),
    "SRC_CF":    ("ClassyFire",          ["InChIKey \\textrightarrow{} parent"]),
    "SRC_CHEBI": ("ChEBI ontology",      ["EBI · OWL"]),
    "SRC_ML":    ("Embeddings + LLM",    ["MiniLM · Claude"]),
    "HUB":       ("Regulatory substances",
                  ["17,523 dbo:ChemicalSubstance / skos:Concept",
                   "CAS · EC number · InChIKey"]),
    "COLL":      ("13 skos:Collection",  ["7 linkability tiers", "6 embedding clusters"]),
    "CHEMONT":   ("ChemOnt 2.1",         ["3,943 skos:Concept",
                                          "4 xkos levels: kingdom \\textrightarrow{} subclass"]),
    "CHEBI_S":   ("ChEBI substances",    ["owl:Class · SMILES · InChI"]),
    "CHEBI_R":   ("ChEBI roles",         ["owl:Class in rdfs:subClassOf", "hierarchy · agrochemical"]),
    "ANNOT":     ("oa:Annotation",       ["48,438 annotations \\ \\ · \\ \\ 40,644 list memberships + 7,794 ChEBI roles"]),
    "SCOPE":     ("Scope mapping",       ["304 SKOS triples", "kept separate"]),
}

FILL = {"src": "volight", "hub": "voyellow", "layer": "white",
        "soft": "volight", "sep": "white"}

# The RDF triples the diagram must assert, as subject -> object.
# An arrow labelled with a predicate always runs from subject to object.
# Verified against substances_taxonomy.ttl (37 MB, 721,603 lines):
#   skos:broader     18,003  csc: -> csc:   (substance -> ChemOnt concept)
#   skos:exactMatch   4,760  csc: -> obo:   (substance -> ChEBI substance class)
#   skos:member          13  collections -> substances
#   oa:hasTarget     48,438  annotation -> substance
#   oa:hasBody       48,438  = 40,644 -> regulatory lists, 7,794 -> ChEBI roles
#   rdfs:subClassOf   5,596  obo: -> obo: ONLY (ChEBI-internal hierarchy)
#   RO_0000087 / owl:Restriction: 0 occurrences
# substances_taxonomy_levels.ttl adds the 4 xkos:ClassificationLevel nodes
# (kingdom/superclass/class/subclass, xkos:depth 1-4) chained as an rdf:List
# off conceptscheme:chemical_substance, each with skos:member to its ChemOnt
# concepts, plus denormalised wk:kingdom/superclass/class/subclass on every
# substance. That file carries NO ChEBI; the taxonomy file carries no levels.
# There is NO substance-class -> role-class edge: the role reaches the
# substance through the annotation layer, not through the ChEBI class tree.
TRIPLES = {
    ("COLL",    "skos:member",     "HUB"),
    ("ANNOT",   "oa:hasTarget",    "HUB"),
    ("ANNOT",   "oa:hasBody",      "CHEBI_R"),
    ("HUB",     "skos:broader",    "CHEMONT"),
    ("HUB",     "skos:exactMatch", "CHEBI_S"),
}


class Layout:
    def __init__(self, w, h, name):
        self.W, self.H, self.name = w, h, name
        self.N = {}          # key -> (x, y, w, h, style)
        self.asserted = set()
        self.E = []          # (pts[list of (x,y)], dashed, label)

    def box(self, key, x, y, w, h, style="layer"):
        self.N[key] = (x, y, w, h, style)

    # ---- anchor-based routing (direction is explicit) ------------------
    def anc(self, key, side, at=None):
        """Point on a box border. `at` fixes the free coordinate absolutely."""
        x, y, w, h = self.N[key][:4]
        if side in "tb":
            xx = x + w/2 if at is None else min(max(at, x+8), x+w-8)
            return (xx, y + h if side == "t" else y)
        yy = y + h/2 if at is None else min(max(at, y+8), y+h-8)
        return (x if side == "l" else x + w, yy)

    def link(self, a, sa, b, sb, mids=(), label="", dashed=False,
             at_a=None, at_b=None):
        """Explicit directed a -> b, orthogonal polyline."""
        if label:
            self.asserted.add((a, label, b))
        # vertical pair with no waypoints: share one x
        if not mids and sa in "tb" and sb in "tb":
            xa = self.anc(a, sa, at_a)[0]
            p = self.anc(a, sa, at_a)
            q = self.anc(b, sb, xa if at_b is None else at_b)
            if abs(p[0]-q[0]) > 1e-6:            # target too narrow: use its x
                p = self.anc(a, sa, q[0])
            self.E.append(([p, (p[0], q[1])] if abs(p[0]-q[0]) < 1e-6
                           else [p, (p[0], (p[1]+q[1])/2), (q[0], (p[1]+q[1])/2), q],
                           dashed, label))
            return
        # horizontal pair with no waypoints: share one y
        if not mids and sa in "lr" and sb in "lr":
            p = self.anc(a, sa, at_a)
            q = self.anc(b, sb, p[1] if at_b is None else at_b)
            if abs(p[1]-q[1]) > 1e-6:
                p = self.anc(a, sa, q[1])
            self.E.append(([p, q], dashed, label))
            return
        p, q = self.anc(a, sa, at_a), self.anc(b, sb, at_b)
        pts = [p] + [tuple(m) for m in mids] + [q]
        out = [pts[0]]
        for nxt in pts[1:]:
            cur = out[-1]
            if abs(cur[0]-nxt[0]) > 1e-6 and abs(cur[1]-nxt[1]) > 1e-6:
                out.append((nxt[0], cur[1]))
            out.append(nxt)
        self.E.append((out, dashed, label))

    # ---- orthogonal routers -------------------------------------------
    def vert(self, a, b, label="", dashed=False, xoff=0):
        """straight vertical: a above b"""
        ax, ay, aw, ah, _ = self.N[a]
        bx, by, bw, bh, _ = self.N[b]
        x = min(max(ax + aw/2 + xoff, bx + 6), bx + bw - 6)
        self.E.append(([(x, ay), (x, by + bh)], dashed, label))

    def horz(self, a, b, label="", dashed=False, yoff=0):
        """straight horizontal: a left of b"""
        ax, ay, aw, ah, _ = self.N[a]
        bx, by, bw, bh, _ = self.N[b]
        y = min(max(ay + ah/2 + yoff, by + 6), by + bh - 6)
        self.E.append(([(ax + aw, y), (bx, y)], dashed, label))

    def elbow_dr(self, a, b, ych, label="", dashed=False):
        """a -> down to channel ych -> across -> down into top of b"""
        ax, ay, aw, ah, _ = self.N[a]
        bx, by, bw, bh, _ = self.N[b]
        x1, x2 = ax + aw/2, bx + bw/2
        self.E.append(([(x1, ay), (x1, ych), (x2, ych), (x2, by + bh)], dashed, label))

    def elbow_rd(self, a, b, xch, label="", dashed=False):
        """a -> right to channel xch -> vertical -> into left of b"""
        ax, ay, aw, ah, _ = self.N[a]
        bx, by, bw, bh, _ = self.N[b]
        y1, y2 = ay + ah/2, by + bh/2
        self.E.append(([(ax + aw, y1), (xch, y1), (xch, y2), (bx, y2)], dashed, label))

    def elbow_under(self, a, b, ych, label="", dashed=False):
        """a -> down below -> across -> up into bottom of b"""
        ax, ay, aw, ah, _ = self.N[a]
        bx, by, bw, bh, _ = self.N[b]
        x1, x2 = ax + aw/2, bx + bw/2
        self.E.append(([(x1, ay), (x1, ych), (x2, ych), (x2, by)], dashed, label))

    # ---- validation ----------------------------------------------------
    def _seg_hits(self, p, q, rect, pad=1.5):
        x, y, w, h = rect
        x0, y0, x1, y1 = x - pad, y - pad, x + w + pad, y + h + pad
        # axis-aligned segments only
        if abs(p[0] - q[0]) < 1e-6:                       # vertical
            xx = p[0]
            if not (x0 < xx < x1): return False
            lo, hi = sorted([p[1], q[1]])
            return not (hi <= y0 or lo >= y1)
        if abs(p[1] - q[1]) < 1e-6:                       # horizontal
            yy = p[1]
            if not (y0 < yy < y1): return False
            lo, hi = sorted([p[0], q[0]])
            return not (hi <= x0 or lo >= x1)
        return False

    def validate(self, verbose=True):
        errs = []
        for k, (x, y, w, h, _) in self.N.items():
            if x < -0.01 or y < -0.01 or x + w > self.W + .01 or y + h > self.H + .01:
                errs.append(f"{k} outside canvas")
        ks = list(self.N)
        for i, a in enumerate(ks):
            for b in ks[i+1:]:
                ax, ay, aw, ah, _ = self.N[a]; bx, by, bw, bh, _ = self.N[b]
                if not (ax+aw <= bx or bx+bw <= ax or ay+ah <= by or by+bh <= ay):
                    errs.append(f"{a} overlaps {b}")
        # edges through boxes
        for idx, (pts, dashed, lab) in enumerate(self.E):
            for j in range(len(pts)-1):
                p, q = pts[j], pts[j+1]
                for k, v in self.N.items():
                    if self._seg_hits(p, q, v[:4]):
                        # allow touching the two endpoints' own boxes
                        if j == 0 and self._is_endpoint(pts[0], v[:4]): continue
                        if j == len(pts)-2 and self._is_endpoint(pts[-1], v[:4]): continue
                        errs.append(f"edge {idx} ({lab or '-'}) passes through {k}")
        wrong = self.asserted - TRIPLES
        for a, pr, b in sorted(wrong):
            rev = (b, pr, a)
            errs.append(f"WRONG DIRECTION {a} -{pr}-> {b}"
                        + ("  (should be reversed)" if rev in TRIPLES else ""))
        for t in sorted(TRIPLES - self.asserted):
            errs.append(f"MISSING {t[0]} -{t[1]}-> {t[2]}")
        if verbose and errs:
            print(f"  !! {self.name}: " + "; ".join(errs[:6]))
        return errs

    def _is_endpoint(self, pt, rect, tol=2.0):
        x, y, w, h = rect
        return (x-tol <= pt[0] <= x+w+tol) and (y-tol <= pt[1] <= y+h+tol)

    def crossings(self):
        segs = []
        for pts, d, l in self.E:
            for j in range(len(pts)-1):
                segs.append((pts[j], pts[j+1]))
        n = 0
        for i in range(len(segs)):
            for j in range(i+1, len(segs)):
                (a, b), (c, d) = segs[i], segs[j]
                av = abs(a[0]-b[0]) < 1e-6
                cv = abs(c[0]-d[0]) < 1e-6
                if av == cv: continue
                if av: vx, vy0, vy1 = a[0], *sorted([a[1], b[1]]); hy, hx0, hx1 = c[1], *sorted([c[0], d[0]])
                else:  vx, vy0, vy1 = c[0], *sorted([c[1], d[1]]); hy, hx0, hx1 = a[1], *sorted([a[0], b[0]])
                if hx0 < vx < hx1 and vy0 < hy < vy1: n += 1
        return n

    # ---- emit ----------------------------------------------------------
    def emit(self, path):
        errs = self.validate()
        if errs:
            raise SystemExit(f"LAYOUT ERRORS in {self.name}:\n  " + "\n  ".join(errs))
        L = ["\\begin{tikzpicture}[x=1mm,y=1mm]",
             f"\\path[use as bounding box] (0,0) rectangle ({self.W},{self.H});"]
        for pts, dashed, lab in self.E:
            style = "dash pattern=on 4mm off 3mm, " if dashed else ""
            path_s = " -- ".join(f"({p[0]:.1f},{p[1]:.1f})" for p in pts)
            L.append(f"\\draw[{style}-{{Stealth[length=5mm,width=3.6mm]}}, line width=1.0mm, "
                     f"vodark, rounded corners=3mm] {path_s};")
            if lab:
                # label where there is most clearance from every box
                def clearance(px, py):
                    m = 1e9
                    for bx, by, bw, bh, _ in self.N.values():
                        dx = max(bx - px, 0, px - (bx + bw))
                        dy = max(by - py, 0, py - (by + bh))
                        m = min(m, (dx*dx + dy*dy) ** 0.5)
                    return m
                best, score = None, -1
                for j in range(len(pts)-1):
                    seglen = abs(pts[j][0]-pts[j+1][0]) + abs(pts[j][1]-pts[j+1][1])
                    if seglen < 6: continue
                    cx = (pts[j][0]+pts[j+1][0])/2
                    cy = (pts[j][1]+pts[j+1][1])/2
                    sc = clearance(cx, cy)
                    if sc > score: score, best = sc, j
                if best is None: best = 0
                mx = (pts[best][0]+pts[best+1][0])/2
                my = (pts[best][1]+pts[best+1][1])/2
                L.append(f"\\node[inner sep=1.0mm, fill=white, text=vodark, "
                         f"font=\\fontsize{{16}}{{19}}\\selectfont\\ttfamily] "
                         f"at ({mx:.1f},{my:.1f}) {{{lab}}};")
        for k, (x, y, w, h, style) in self.N.items():
            title, det = TEXT[k]
            dash = ", dash pattern=on 5mm off 3.5mm" if style == "sep" else ""
            body = f"\\fontsize{{25}}{{28}}\\selectfont\\bfseries {title}"
            if det:
                body += ("\\\\[1.2mm]\\fontsize{18}{22}\\selectfont\\mdseries "
                         + "\\\\ ".join(det))
            L.append(f"\\node[anchor=south west, minimum width={w}mm, minimum height={h}mm, "
                     f"text width={w-10}mm, align=center, inner sep=0pt, fill={FILL[style]}, "
                     f"draw=vodark, line width=0.8mm{dash}, "
                     f"font=\\fontsize{{25}}{{28}}\\selectfont] at ({x},{y}) {{{body}}};")
        L.append("\\end{tikzpicture}%")
        open(path, "w").write("\n".join(L) + "\n")
        return len(self.N), len(self.E), self.crossings()


# ==================================================================== OPTION 1
def option1():
    """TIERS: core band across the full width; sources above, facets below.
    A left-hand channel carries the annotation relation back up to the core."""
    L = Layout(780, 322, "option1-tiers")
    CH = 22.0                                   # left return channel
    src_x = [CH, 176, 330, 484, 638]
    src_w = [146, 146, 146, 146, 142]
    for k, x, w in zip(["SRC_REG", "SRC_PUB", "SRC_CF", "SRC_CHEBI", "SRC_ML"],
                       src_x, src_w):
        L.box(k, x, 258, w, 64, "src")
    L.box("HUB", CH, 152, 758, 72, "hub")
    L.box("COLL",    CH,  46, 140, 72, "layer")
    L.box("CHEMONT", 176, 46, 140, 72, "layer")
    L.box("CHEBI_S", 330, 46, 140, 72, "layer")
    L.box("CHEBI_R", 520, 46, 120, 72, "layer")
    L.box("SCOPE",   654, 46, 126, 72, "sep")
    L.box("ANNOT",   CH,   0, 758, 34, "soft")

    for k, x, w in zip(["SRC_REG", "SRC_PUB", "SRC_CF", "SRC_CHEBI", "SRC_ML"],
                       src_x, src_w):
        L.link(k, "b", "HUB", "t", at_a=x+w/2, at_b=x+w/2)
    L.link("HUB", "b", "CHEMONT", "t", at_a=246, at_b=246, label="skos:broader")
    L.link("HUB", "b", "CHEBI_S", "t", at_a=400, at_b=400, label="skos:exactMatch")
    L.link("COLL", "t", "HUB", "b", at_a=92, at_b=92, label="skos:member")
    L.link("HUB", "b", "SCOPE", "t", at_a=717, at_b=717, dashed=True)
    L.link("ANNOT", "t", "CHEBI_R", "b", at_a=580, at_b=580, label="oa:hasBody")
    L.link("ANNOT", "l", "HUB", "l", mids=[(11, 17), (11, 188)],
           label="oa:hasTarget", at_a=17, at_b=188)
    return L


# ==================================================================== OPTION 2
def option2():
    """FAN: core in the middle. Incoming relations (member / hasTarget) sit
    above and below the core; outgoing ones fan out to the right."""
    L = Layout(780, 322, "option2-fan")
    entry = [128, 149, 170, 191, 206]
    for i, (k, ey) in enumerate(zip(
            ["SRC_REG", "SRC_PUB", "SRC_CF", "SRC_CHEBI", "SRC_ML"], entry)):
        L.box(k, 0, 258 - i*62, 150, 52, "src")
    L.box("HUB",   196, 118, 190, 96, "hub")
    L.box("COLL",  216, 240, 150, 68, "layer")
    L.box("ANNOT", 216,  20, 150, 68, "soft")
    L.box("CHEMONT", 448, 240, 150, 68, "layer")
    L.box("CHEBI_S", 448,  80, 150, 68, "layer")
    L.box("CHEBI_R", 650,  80, 130, 68, "layer")
    L.box("SCOPE",   650, 240, 130, 68, "sep")

    for k, ey in zip(["SRC_REG", "SRC_PUB", "SRC_CF", "SRC_CHEBI", "SRC_ML"], entry):
        L.link(k, "r", "HUB", "l", mids=[(173, ey)], at_b=ey)
    L.link("COLL",  "b", "HUB", "t", at_a=291, at_b=291, label="skos:member")
    L.link("ANNOT", "t", "HUB", "b", at_a=291, at_b=291, label="oa:hasTarget")
    L.link("HUB", "r", "CHEMONT", "l", mids=[(415, 274)], at_a=190, label="skos:broader")
    L.link("HUB", "r", "CHEBI_S", "l", mids=[(415, 114)], at_a=142, label="skos:exactMatch")
    L.link("ANNOT", "r", "CHEBI_R", "b", mids=[(715, 54)], at_a=54, at_b=715,
           label="oa:hasBody")
    L.link("HUB", "t", "SCOPE", "t", mids=[(210, 315), (715, 315)],
           at_a=210, at_b=715, dashed=True)
    return L


# ==================================================================== OPTION 3
def option3():
    """PIPELINE: strict left-to-right stages; the two inverse relations
    (member / hasTarget) return through their own channel."""
    L = Layout(780, 340, "option3-pipeline")
    for i, k in enumerate(["SRC_REG", "SRC_PUB", "SRC_CF", "SRC_CHEBI", "SRC_ML"]):
        L.box(k, 0, 282 - i*60, 148, 50, "src")
    L.box("HUB", 214, 130, 172, 110, "hub")
    L.box("CHEMONT", 448, 244, 150, 68, "layer")
    L.box("CHEBI_S", 448, 166, 150, 68, "layer")
    L.box("COLL",    448,  88, 150, 68, "layer")
    L.box("ANNOT",   448,  10, 150, 58, "soft")
    L.box("CHEBI_R", 650, 166, 130, 68, "layer")
    L.box("SCOPE",   650, 244, 130, 68, "sep")

    for k in ["SRC_REG", "SRC_PUB", "SRC_CF", "SRC_CHEBI", "SRC_ML"]:
        L.link(k, "r", "HUB", "l", mids=[(186, 0)])
    L.link("HUB", "r", "CHEMONT", "l", mids=[(410, 278)], label="skos:broader", at_a=217)
    L.link("HUB", "r", "CHEBI_S", "l", mids=[(410, 200)], label="skos:exactMatch", at_a=183)
    L.link("COLL",  "l", "HUB", "r", mids=[(428, 155)], label="skos:member", at_b=155)
    L.link("ANNOT", "l", "HUB", "r", mids=[(420, 140)], label="oa:hasTarget", at_b=140)
    L.link("ANNOT", "r", "CHEBI_R", "b", mids=[(715, 39)], label="oa:hasBody", at_a=39)
    L.link("HUB", "t", "SCOPE", "t", mids=[(300, 326), (715, 326)], dashed=True, at_a=300, at_b=715)
    return L


if __name__ == "__main__":
    for fn, out in [(option1, "kg_opt1.tex"), (option2, "kg_opt2.tex"),
                    (option3, "kg_opt3.tex")]:
        L = fn()
        n, e, x = L.emit(out)
        print(f"{out:14s} {L.name:20s} {n:2d} nodes, {e:2d} edges, "
              f"{x} crossings, canvas {L.W}x{L.H} mm  -> VALID")

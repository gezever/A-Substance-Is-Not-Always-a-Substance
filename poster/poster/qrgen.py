#!/usr/bin/env python3
"""Minimal QR encoder (byte mode, version 5, EC level M) with self-verification.
Version 5-M: 37x37 modules, 2 blocks x (43 data + 24 EC) codewords = 134 total."""

URL = "https://github.com/gezever/A-Substance-Is-Not-Always-a-Substance"

VER, SIZE = 5, 37
NBLOCKS, DATA_PER_BLOCK, EC_PER_BLOCK = 2, 43, 24
TOTAL_DATA = NBLOCKS * DATA_PER_BLOCK          # 86
ALIGN_CENTERS = [6, 30]
ECL_BITS = 0b00                                # level M
REMAINDER_BITS = 7

# ----------------------------------------------------------------- GF(256)
EXP = [0]*512
LOG = [0]*256
x = 1
for i in range(255):
    EXP[i] = x
    LOG[x] = i
    x <<= 1
    if x & 0x100:
        x ^= 0x11D
for i in range(255, 512):
    EXP[i] = EXP[i-255]

def gmul(a, b):
    if a == 0 or b == 0:
        return 0
    return EXP[LOG[a] + LOG[b]]

def rs_generator(n):
    g = [1]
    for i in range(n):
        g2 = [0]*(len(g)+1)
        for j, c in enumerate(g):
            g2[j] ^= c
            g2[j+1] ^= gmul(c, EXP[i])
        g = g2
    return g

def rs_encode(data, n):
    gen = rs_generator(n)
    res = list(data) + [0]*n
    for i in range(len(data)):
        coef = res[i]
        if coef:
            for j in range(1, len(gen)):
                res[i+j] ^= gmul(gen[j], coef)
    return res[len(data):]

def rs_syndromes(codeword, n):
    """All syndromes must be 0 for a valid RS codeword."""
    out = []
    for i in range(n):
        s = 0
        for c in codeword:
            s = gmul(s, EXP[i]) ^ c
        out.append(s)
    return out

# ------------------------------------------------------------ data encoding
def encode_data(text):
    b = text.encode("utf-8")
    assert len(b) < 256
    bits = "0100" + format(len(b), "08b") + "".join(format(c, "08b") for c in b)
    cap = TOTAL_DATA * 8
    assert len(bits) <= cap, f"payload too long: {len(bits)} > {cap}"
    bits += "0" * min(4, cap - len(bits))                 # terminator
    if len(bits) % 8:
        bits += "0" * (8 - len(bits) % 8)                 # byte align
    pad = [0xEC, 0x11]
    i = 0
    while len(bits) < cap:
        bits += format(pad[i % 2], "08b")
        i += 1
    return [int(bits[i:i+8], 2) for i in range(0, len(bits), 8)]

def build_codewords(text):
    data = encode_data(text)
    blocks = [data[i*DATA_PER_BLOCK:(i+1)*DATA_PER_BLOCK] for i in range(NBLOCKS)]
    ecs = [rs_encode(b, EC_PER_BLOCK) for b in blocks]
    # verify every block is a genuine RS codeword
    for b, e in zip(blocks, ecs):
        assert all(s == 0 for s in rs_syndromes(b + e, EC_PER_BLOCK)), "RS syndrome != 0"
    out = []
    for i in range(DATA_PER_BLOCK):
        for b in blocks:
            out.append(b[i])
    for i in range(EC_PER_BLOCK):
        for e in ecs:
            out.append(e[i])
    return out, blocks

# --------------------------------------------------------------- matrix
def new_matrix():
    return [[None]*SIZE for _ in range(SIZE)]

def place_finder(m, r, c):
    for dr in range(-1, 8):
        for dc in range(-1, 8):
            rr, cc = r+dr, c+dc
            if not (0 <= rr < SIZE and 0 <= cc < SIZE):
                continue
            if dr in (-1, 7) or dc in (-1, 7):
                m[rr][cc] = 0
            elif dr in (0, 6) or dc in (0, 6) or (2 <= dr <= 4 and 2 <= dc <= 4):
                m[rr][cc] = 1
            else:
                m[rr][cc] = 0

def place_function_patterns(m):
    place_finder(m, 0, 0)
    place_finder(m, 0, SIZE-7)
    place_finder(m, SIZE-7, 0)
    # timing
    for i in range(SIZE):
        if m[6][i] is None:
            m[6][i] = 1 if i % 2 == 0 else 0
        if m[i][6] is None:
            m[i][6] = 1 if i % 2 == 0 else 0
    # alignment
    for r in ALIGN_CENTERS:
        for c in ALIGN_CENTERS:
            if (r, c) in ((6, 6), (6, ALIGN_CENTERS[-1]), (ALIGN_CENTERS[-1], 6)):
                continue
            for dr in range(-2, 3):
                for dc in range(-2, 3):
                    m[r+dr][c+dc] = 1 if (abs(dr) == 2 or abs(dc) == 2 or (dr == 0 and dc == 0)) else 0
    # dark module
    m[4*VER+9][8] = 1

def format_positions():
    """The two copies of the 15-bit format information."""
    a = [(8, i) for i in range(6)] + [(8, 7), (8, 8), (7, 8)] + [(i, 8) for i in range(5, -1, -1)]
    b = [(SIZE-1-i, 8) for i in range(7)] + [(8, SIZE-8+i) for i in range(8)]
    return a, b

def reserve_format(m):
    for pos in format_positions():
        for r, c in pos:
            m[r][c] = 0

def data_positions(m):
    """Zigzag order of free modules, right to left, skipping the vertical timing column."""
    pos = []
    col = SIZE - 1
    upward = True
    while col > 0:
        if col == 6:
            col -= 1
        rows = range(SIZE-1, -1, -1) if upward else range(SIZE)
        for r in rows:
            for c in (col, col-1):
                if m[r][c] is None:
                    pos.append((r, c))
        upward = not upward
        col -= 2
    return pos

MASKS = [
    lambda r, c: (r + c) % 2 == 0,
    lambda r, c: r % 2 == 0,
    lambda r, c: c % 3 == 0,
    lambda r, c: (r + c) % 3 == 0,
    lambda r, c: (r//2 + c//3) % 2 == 0,
    lambda r, c: (r*c) % 2 + (r*c) % 3 == 0,
    lambda r, c: ((r*c) % 2 + (r*c) % 3) % 2 == 0,
    lambda r, c: ((r+c) % 2 + (r*c) % 3) % 2 == 0,
]

def penalty(m):
    score = 0
    # rule 1: runs of 5+
    for line in list(m) + [list(col) for col in zip(*m)]:
        run, prev = 1, line[0]
        for v in line[1:]:
            if v == prev:
                run += 1
            else:
                if run >= 5:
                    score += 3 + (run - 5)
                run, prev = 1, v
        if run >= 5:
            score += 3 + (run - 5)
    # rule 2: 2x2 blocks
    for r in range(SIZE-1):
        for c in range(SIZE-1):
            if m[r][c] == m[r][c+1] == m[r+1][c] == m[r+1][c+1]:
                score += 3
    # rule 3: finder-like patterns
    pat1 = [1,0,1,1,1,0,1,0,0,0,0]
    pat2 = [0,0,0,0,1,0,1,1,1,0,1]
    for line in list(m) + [list(col) for col in zip(*m)]:
        for i in range(SIZE-10):
            seg = line[i:i+11]
            if seg == pat1 or seg == pat2:
                score += 40
    # rule 4: dark ratio
    dark = sum(sum(row) for row in m)
    pct = dark * 100 / (SIZE*SIZE)
    score += 10 * (int(abs(pct-50)//5))
    return score

def format_bits(ecl, mask):
    v = (ecl << 3) | mask
    d = v << 10
    g = 0b10100110111
    for i in range(4, -1, -1):
        if d & (1 << (i+10)):
            d ^= g << i
    return ((v << 10) | d) ^ 0b101010000010010

def build(text):
    codewords, blocks = build_codewords(text)
    bits = "".join(format(c, "08b") for c in codewords) + "0"*REMAINDER_BITS

    base = new_matrix()
    place_function_patterns(base)
    reserve_format(base)
    positions = data_positions(base)
    assert len(positions) == len(bits), f"{len(positions)} slots vs {len(bits)} bits"

    best = None
    for mask in range(8):
        m = [row[:] for row in base]
        for (r, c), bit in zip(positions, bits):
            m[r][c] = int(bit) ^ (1 if MASKS[mask](r, c) else 0)
        fb = format_bits(ECL_BITS, mask)
        fa, fbpos = format_positions()
        for i, (r, c) in enumerate(fa):
            m[r][c] = (fb >> (14-i)) & 1
        for i, (r, c) in enumerate(fbpos):
            m[r][c] = (fb >> (14-i)) & 1
        m[4*VER+9][8] = 1
        p = penalty(m)
        if best is None or p < best[0]:
            best = (p, mask, m)
    return best[2], best[1], best[0], positions, codewords

# --------------------------------------------------------------- verify
def read_back(m, mask, positions, expected_text):
    bits = ""
    for (r, c) in positions:
        bits += str(m[r][c] ^ (1 if MASKS[mask](r, c) else 0))
    cw = [int(bits[i:i+8], 2) for i in range(0, len(bits)//8*8, 8)][:NBLOCKS*(DATA_PER_BLOCK+EC_PER_BLOCK)]
    # de-interleave data codewords
    blocks = [[] for _ in range(NBLOCKS)]
    idx = 0
    for i in range(DATA_PER_BLOCK):
        for b in range(NBLOCKS):
            blocks[b].append(cw[idx]); idx += 1
    data = sum(blocks, [])
    dbits = "".join(format(c, "08b") for c in data)
    assert dbits[:4] == "0100", "mode indicator wrong"
    n = int(dbits[4:12], 2)
    payload = bytes(int(dbits[12+8*i:20+8*i], 2) for i in range(n))
    got = payload.decode("utf-8")
    assert got == expected_text, f"round-trip mismatch:\n  {got!r}\n  {expected_text!r}"
    return got

if __name__ == "__main__":
    m, mask, pen, positions, cw = build(URL)
    got = read_back(m, mask, positions, URL)
    print("payload      :", got)
    print("version/ecl  : %d-M   size %dx%d" % (VER, SIZE, SIZE))
    print("mask chosen  :", mask, " penalty", pen)
    print("RS syndromes : all zero (checked per block)")
    print("round-trip   : OK")

    # quiet zone 4 modules, emit TikZ (vector) and PNG (preview)
    Q = 4
    with open("qrcode.tex", "w") as f:
        f.write("%% QR code for %s -- generated, version %d-M\n" % (URL, VER))
        f.write("\\begin{tikzpicture}[x=1mm,y=1mm]\n")
        f.write("\\path[use as bounding box] (0,0) rectangle (%d,%d);\n" % (SIZE+2*Q, SIZE+2*Q))
        f.write("\\fill[white] (0,0) rectangle (%d,%d);\n" % (SIZE+2*Q, SIZE+2*Q))
        f.write("\\begin{scope}[shift={(%d,%d)}]\n" % (Q, Q))
        for r in range(SIZE):
            for c in range(SIZE):
                if m[r][c]:
                    f.write("\\fill (%d,%d) rectangle ++(1,1);\n" % (c, SIZE-1-r))
        f.write("\\end{scope}\n\\end{tikzpicture}%\n")

    from PIL import Image
    S = 12
    img = Image.new("RGB", ((SIZE+2*Q)*S, (SIZE+2*Q)*S), (255, 255, 255))
    px = img.load()
    for r in range(SIZE):
        for c in range(SIZE):
            if m[r][c]:
                for dy in range(S):
                    for dx in range(S):
                        px[(c+Q)*S+dx, (r+Q)*S+dy] = (0, 0, 0)
    img.save("qrcode.png")
    print("written      : qrcode.tex (vector), qrcode.png")

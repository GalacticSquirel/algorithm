"""The reduction in README.md ("Finding the minimum is NP-hard"), as a program.

Input: a planar graph of maximum degree 3 with an orthogonal grid drawing (vertices
at grid points, edges as rectilinear polylines). Output: an emitter instance in the
assignment's input format, together with the number of forced cells and of
subdivisions, so that

    min emitters = forced + VC(G) + subdivisions.

Cells (i,j) are the odd nodes (2i+1, 2j+1); corners (i,j) are the even nodes
(2i, 2j). Every even node holds a parasite. The asserts check the three gadget
invariants (I1)-(I3) of the README on every generated instance.

Run: python3 gen.py   (writes k4.txt, c5.txt, star.txt and prints the expected minimum)
Not part of the Java submission."""
import itertools, sys

S = 32     # scale of the coarse drawing (enough room for every gadget)
MARGIN = 4

def sgn(a): return (a > 0) - (a < 0)

def expand(poly):
    cells = [poly[0]]
    for (x0, y0), (x1, y1) in zip(poly, poly[1:]):
        dx, dy = sgn(x1 - x0), sgn(y1 - y0)
        assert dx == 0 or dy == 0
        x, y = x0, y0
        while (x, y) != (x1, y1):
            x, y = x + dx, y + dy
            cells.append((x, y))
    return cells

def build(vertices, edges):
    # scale
    P = {v: (S * x + MARGIN, S * y + MARGIN) for v, (x, y) in vertices.items()}
    chains = []
    for (u, v, poly) in edges:
        pts = [(S * x + MARGIN, S * y + MARGIN) for x, y in poly]
        assert pts[0] == P[u] and pts[-1] == P[v]
        chains.append([u, v, expand(pts)])
    # degree-3 gadget: T-shape with arms w, e collinear and n perpendicular
    arms = {v: [] for v in vertices}
    for idx, (u, v, c) in enumerate(chains):
        arms[u].append((idx, 0, (c[1][0] - c[0][0], c[1][1] - c[0][1])))
        arms[v].append((idx, 1, (c[-2][0] - c[-1][0], c[-2][1] - c[-1][1])))
    for v, a in arms.items():
        assert len(a) <= 3, v
        if len(a) < 3:
            continue
        dirs = [d for _, _, d in a]
        # find collinear pair
        for i, j in itertools.combinations(range(3), 2):
            if dirs[i][0] == -dirs[j][0] and dirs[i][1] == -dirs[j][1]:
                k = 3 - i - j
                break
        else:
            raise ValueError("vertex %s is not a T" % v)
        e = dirs[i]; n = dirs[k]; s = (-n[0], -n[1])
        idx, end, _ = a[i]
        c = chains[idx][2]
        if end == 1:
            c.reverse()
        x, y = c[0]
        # c[1], c[2] are x+e, x+2e; replace by x+e+s, x+2e+s, x+3e+s, then x+3e (kept)
        for t in (1, 2, 3, 4):
            assert c[t] == (x + t * e[0], y + t * e[1]), "arm too short"
        # diagonal step away from n, two straight cells, diagonal step back
        new = [c[0]] + [(x + t * e[0] + s[0], y + t * e[1] + s[1]) for t in (1, 2, 3)] + c[4:]
        if end == 1:
            new.reverse()
        chains[idx][2] = new
    # parity: internal cells must be even; fix with a diagonal step + jog on a straight run
    for ch in chains:
        c = ch[2]
        if (len(c) - 2) % 2 == 1:
            ch[2] = parity_fix(c)
        assert (len(ch[2]) - 2) % 2 == 0
    return P, chains

def parity_fix(c):
    # On a straight run along d: three diagonal steps towards n, two straight
    # cells, then straight back across three rows. Three diagonal steps flip
    # the parity of the cell count.
    L = 14
    for start in range(3, len(c) - L - 4):
        run = c[start - 3:start + L + 3]
        d = (run[1][0] - run[0][0], run[1][1] - run[0][1])
        if all((run[t + 1][0] - run[t][0], run[t + 1][1] - run[t][1]) == d for t in range(len(run) - 1)):
            run = c[start:]
            n = (-d[1], d[0])
            x, y = run[0]
            at = lambda a, b: (x + a * d[0] + b * n[0], y + a * d[1] + b * n[1])
            detour = [at(1, 1), at(2, 2), at(3, 3), at(4, 3), at(5, 3), at(6, 3),
                      at(6, 2), at(6, 1), at(6, 0)]
            return c[:start + 1] + detour + c[start + 7:]
    raise ValueError("no straight run for parity fix")

def instance(vertices, edges):
    P, chains = build(vertices, edges)
    V = set()
    owner = {}
    for idx, (_, _, c) in enumerate(chains):
        for t, cell in enumerate(c):
            ends = t in (0, len(c) - 1)
            if cell in owner:
                assert ends and owner[cell] == "end", ("chains overlap", cell)
            owner[cell] = "end" if ends else idx
            V.add(cell)
    W = max(max(x for x, _ in V), max(y for _, y in V)) + MARGIN + 1  # cells 0..W-1
    l = W
    var = lambda c: c in V
    def block(p):  # cells around corner p (inside board)
        i, j = p
        return [(a, b) for a in (i - 1, i) for b in (j - 1, j) if 0 <= a < l and 0 <= b < l]
    pairs = []
    for _, _, c in chains:
        pairs += list(zip(c, c[1:]))
    # intended corner for each pair
    intended = {}
    for a, b in pairs:
        near = [(a[0] + u, a[1] + w) for u in (0, 1) for w in (0, 1)]
        common = [p for p in near if a in block(p) and b in block(p)]
        good = [p for p in common if set(q for q in block(p) if var(q)) == {a, b}]
        assert good, ("pair has no clean corner", a, b)
        intended[good[0]] = (a, b)
    demand = {}
    for p in itertools.product(range(l + 1), repeat=2):
        f = sum(1 for q in block(p) if not var(q))
        assert f >= 1, ("corner with only emitter-choice cells", p)
        demand[p] = f + (1 if p in intended else 0)
    # (ii) every non-var cell is in an all-non-var block
    for c in itertools.product(range(l), repeat=2):
        if c in V:
            continue
        i, j = c
        corners = [(i, j), (i + 1, j), (i, j + 1), (i + 1, j + 1)]
        assert any(all(not var(q) for q in block(p)) for p in corners), ("cell not forced", c)
    F = l * l - len(V)
    s = sum((len(c) - 2) // 2 for _, _, c in chains)
    return l, demand, F, s, V, pairs

def write(fn, l, demand):
    with open(fn, "w") as out:
        out.write(f"{l}\n")
        for (i, j), d in sorted(demand.items()):
            out.write(f"{2*i} {2*j} {d}\n")

def vc(vertices, edges):
    vs = list(vertices)
    for k in range(len(vs) + 1):
        for C in itertools.combinations(vs, k):
            if all(u in C or v in C for u, v, _ in edges):
                return k

GRAPHS = {
    # K4: outer A,B,C with D inside
    "k4": ({"A": (2, 0), "B": (4, 2), "C": (2, 4), "D": (2, 2)},
           [("A", "D", [(2, 0), (2, 2)]), ("D", "C", [(2, 2), (2, 4)]), ("D", "B", [(2, 2), (4, 2)]),
            ("A", "C", [(2, 0), (0, 0), (0, 4), (2, 4)]), ("A", "B", [(2, 0), (4, 0), (4, 2)]),
            ("B", "C", [(4, 2), (4, 4), (2, 4)])]),
    # a path a-b-c and a triangle-free cycle of length 5
    "c5": ({"a": (0, 0), "b": (3, 0), "c": (3, 3), "d": (1, 5), "e": (0, 3)},
           [("a", "b", [(0, 0), (3, 0)]), ("b", "c", [(3, 0), (3, 3)]), ("c", "d", [(3, 3), (3, 5), (1, 5)]),
            ("d", "e", [(1, 5), (0, 5), (0, 3)]), ("e", "a", [(0, 3), (0, 0)])]),
    # star K_{1,3} plus a leaf path
    "star": ({"o": (2, 2), "p": (0, 2), "q": (4, 2), "r": (2, 4), "t": (4, 4)},
             [("o", "p", [(2, 2), (0, 2)]), ("o", "q", [(2, 2), (4, 2)]), ("o", "r", [(2, 2), (2, 4)]),
              ("q", "t", [(4, 2), (4, 4)])]),
}

if __name__ == "__main__":
    for name, (vs, es) in GRAPHS.items():
        l, demand, F, s, V, pairs = instance(vs, es)
        write(f"{name}.txt", l, demand)
        k = vc(vs, es)
        print(name, "l", l, "forced", F, "subdivisions", s, "VC(G)", k, "expected min", F + k + s)

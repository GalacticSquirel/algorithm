"""Independent end-to-end check of the reduction, using only the brief's rules.

check.py trusts Lemmas 1 and 2 (it treats the instance as a cell multicover problem).
This script does not. It decides whether a set of emitters works the way the brief
defines it: vertex-disjoint pulses, one pulse per incident edge of an emitter, no
pulse through another emitter or parasite. It does this with its own node-split max
flow (Dinitz), written separately from the Java solver.

1. Lemma 1: on random small grids where every even node is a parasite, max flow
   agrees with "every corner has >= d chosen cells around it".
2. For each graph, at a reduced scale (S = 14, so it runs in minutes):
   - a cell is *necessary* if removing just that cell from "every cell is an emitter"
     makes max flow fall short. Necessary cells are in every solution.
   - a pair {a, b} is a *pair constraint* if removing both makes max flow fall short.
     Every solution contains a or b.
   Then it checks that the pair constraints form exactly G with every edge replaced
   by a path with an even number of inner vertices, computes the minimum vertex cover
   of that graph by DP (not via Lemma 3), and gives the lower bound
   necessary + VC. It builds an emitter set of that size and checks with max flow
   that it works. So the minimum is proven from the rules alone, and it is compared
   with gen.py's F + VC(G) + subdivisions.

Run: python3 semantic.py          (about 15-30 minutes on 4 cores)
     python3 semantic.py K23      (one graph)
Not part of the Java submission."""
import sys, itertools, random
from collections import deque
from multiprocessing import Pool
import gen

S_SMALL = 14


# ---------- the grid, as the brief defines it ----------
def nodes(l):
    return [(x, y) for x in range(2 * l + 1) for y in range(2 * l + 1) if x % 2 == y % 2]


def nbrs(v, l):
    x, y = v
    steps = [(1, 1), (1, -1), (-1, 1), (-1, -1)]
    if x % 2 == 0:                       # even nodes also have the grid-line edges
        steps += [(2, 0), (-2, 0), (0, 2), (0, -2)]
    return [(x + dx, y + dy) for dx, dy in steps if 0 <= x + dx <= 2 * l and 0 <= y + dy <= 2 * l]


class Dinic:
    def __init__(self, n):
        self.g = [[] for _ in range(n)]

    def add(self, u, v, c):
        self.g[u].append([v, c, len(self.g[v])])
        self.g[v].append([u, 0, len(self.g[u]) - 1])

    def flow(self, s, t):
        total = 0
        while True:
            level = [-1] * len(self.g)
            level[s] = 0
            q = deque([s])
            while q:
                u = q.popleft()
                for v, c, _ in self.g[u]:
                    if c > 0 and level[v] < 0:
                        level[v] = level[u] + 1
                        q.append(v)
            if level[t] < 0:
                return total
            it = [0] * len(self.g)

            def dfs(u, f):
                if u == t:
                    return f
                while it[u] < len(self.g[u]):
                    e = self.g[u][it[u]]
                    v, c, r = e
                    if c > 0 and level[v] == level[u] + 1:
                        d = dfs(v, min(f, c))
                        if d:
                            e[1] -= d
                            self.g[v][r][1] += d
                            return d
                    it[u] += 1
                return 0

            while True:
                f = dfs(s, 1 << 30)
                if not f:
                    break
                total += f


def works(l, dem, emitters):
    """True if this emitter set neutralises every parasite (brief, section 3.1)."""
    V = nodes(l)
    idx = {v: i for i, v in enumerate(V)}
    n = len(V)
    D = Dinic(2 * n + 2)
    src, sink = 2 * n, 2 * n + 1
    for v in V:
        i = idx[v]
        if v in dem:                          # parasite: pulses end here
            D.add(i, sink, dem[v])
            continue
        if v in emitters:                     # emitter: one pulse per incident edge
            D.add(src, n + i, len(nbrs(v, l)))
        else:                                 # free node: at most one pulse passes
            D.add(i, n + i, 1)
        for w in nbrs(v, l):
            if w not in emitters:             # never into another emitter
                D.add(n + i, idx[w], 1)
    return D.flow(src, sink) == sum(dem.values())


def lemma1(trials=300):
    rnd = random.Random(1)
    for _ in range(trials):
        l = rnd.randint(1, 4)
        dem = {}
        for v in nodes(l):
            if v[0] % 2 == 0:
                dem[v] = rnd.randint(1, sum(1 for w in nbrs(v, l) if w[0] % 2))
        emit = {v for v in nodes(l) if v[0] % 2 and rnd.random() < 0.8}
        cover = all(sum(1 for w in nbrs(p, l) if w in emit) >= d for p, d in dem.items())
        assert works(l, dem, emit) == cover, (l, dem, emit)
    print("Lemma 1: max flow agrees with the corner condition on", trials, "random grids")


# ---------- end-to-end check of one generated instance ----------
def _necessary(args):
    l, dem, cells, chunk = args
    return [c for c in chunk if not works(l, dem, cells - {c})]


def path_from(adj, c, w, ends):
    p = [c, w]
    while p[-1] not in ends:
        p.append(next(x for x in adj[p[-1]] if x != p[-2]))
    return p


def check(name, vs, es, pool):
    gen.S = S_SMALL
    l, demand, F, s, _, _ = gen.instance(vs, es)
    dem = {(2 * i, 2 * j): d for (i, j), d in demand.items()}
    cells = {v for v in nodes(l) if v[0] % 2}
    assert works(l, dem, cells), "instance must be solvable"
    order = sorted(cells)
    chunks = [(l, dem, cells, order[i::32]) for i in range(32)]
    necessary = set(sum(pool.map(_necessary, chunks), []))
    rest = cells - necessary
    cand = set()
    for p in dem:
        around = sorted(w for w in nbrs(p, l) if w in rest)
        cand |= set(itertools.combinations(around, 2))
    cand = sorted(cand)
    bad = pool.map(_necessary_pair, [(l, dem, cells, cand[i::32]) for i in range(32)])
    pairs = set(sum(bad, []))
    adj = {}
    for a, b in pairs:
        adj.setdefault(a, set()).add(b)
        adj.setdefault(b, set()).add(a)

    # the pair constraints are G with every edge subdivided an even number of times
    vcell = {v: (2 * (gen.S * x + gen.MARGIN) + 1, 2 * (gen.S * y + gen.MARGIN) + 1) for v, (x, y) in vs.items()}
    ends = set(vcell.values())
    cell2v = {c: v for v, c in vcell.items()}
    assert all(len(adj[c]) == 2 for c in adj if c not in ends)
    traced = []
    for v, c in vcell.items():
        assert len(adj.get(c, ())) == sum((e[0] == v) + (e[1] == v) for e in es), v
        for w in adj.get(c, ()):
            p = path_from(adj, c, w, ends)
            assert (len(p) - 2) % 2 == 0
            if c < p[-1]:
                traced.append((tuple(sorted((v, cell2v[p[-1]]))), p))
    assert sorted(t for t, _ in traced) == sorted(tuple(sorted(e[:2])) for e in es)

    # exact minimum vertex cover of the pair graph: try every choice on the vertex
    # cells, then each path is solved by DP
    best, witness = None, None
    vlist = sorted(ends)
    for bits in itertools.product((0, 1), repeat=len(vlist)):
        chosen = {c for c, b in zip(vlist, bits) if b}
        total, pick, ok = len(chosen), set(chosen), True
        for _, p in traced:
            a, b, inner = p[0] in chosen, p[-1] in chosen, p[1:-1]
            # dp[x] = (cost, picks) with x = whether the last vertex is chosen
            dp = {1: (0, [])} if a else {0: (0, [])}
            for cell in inner:
                nd = {}
                for x in (0, 1):
                    opts = [(cv + x, pk + ([cell] if x else [])) for px, (cv, pk) in dp.items() if px or x]
                    if opts:
                        nd[x] = min(opts, key=lambda o: o[0])
                dp = nd
            opts = [val for px, val in dp.items() if px or b]
            if not opts:
                ok = False
                break
            cost, pk = min(opts, key=lambda o: o[0])
            total += cost
            pick |= set(pk)
        if ok and (best is None or total < best):
            best, witness = total, pick
    lower = len(necessary) + best
    emitters = necessary | witness
    assert len(emitters) == lower and works(l, dem, emitters)
    expect = F + vc(vs, es) + s
    print(f"{name}: l={l}, minimum by max flow = {lower}, gen.py predicts {expect}")
    assert lower == expect


def _necessary_pair(args):
    l, dem, cells, chunk = args
    return [ab for ab in chunk if not works(l, dem, cells - set(ab))]


def vc(vs, es):
    for k in range(len(vs) + 1):
        for C in itertools.combinations(vs, k):
            if all(u in C or v in C for u, v, _ in es):
                return k


if __name__ == "__main__":
    lemma1()
    graphs = {k: v for k, v in gen.GRAPHS.items() if k != "Q3"}   # Q3 needs scale >= 27
    if len(sys.argv) > 1:
        graphs = {k: gen.GRAPHS[k] for k in sys.argv[1:]}
    with Pool() as pool:
        for name, (vs, es) in graphs.items():
            check(name, vs, es, pool)

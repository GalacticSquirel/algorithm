"""Independent check of a generated instance, from the file alone.
1. Every odd node is free and every even node a parasite, so pulses are single edges and the
   problem is: choose cells so that each corner p has >= d(p) chosen cells around it.
2. Corners with d = number of cells force those cells.
3. The residual constraints must all be 'at least 1 of these 2 cells'; compute the vertex cover
   of that graph exactly (branch on vertices of degree >= 3, paths by DP).
Run: python3 check.py k4.txt   (prints the minimum, to compare with gen.py)
Not part of the Java submission."""
import sys, itertools

def load(fn):
    nums = list(map(int, open(fn).read().split()))
    l = nums[0]; d = {}
    for i in range(1, len(nums), 3):
        d[(nums[i], nums[i + 1])] = nums[i + 2]
    return l, d

def cells_of(p, l):
    x, y = p
    return [(x + a, y + b) for a in (-1, 1) for b in (-1, 1) if 0 < x + a < 2 * l and 0 < y + b < 2 * l]

def vc_exact(adj):
    # vertices of degree >= 3 branched; the rest is paths/cycles, solved by brute DP per component
    big = [v for v in adj if len(adj[v]) >= 3]
    best = None
    for bits in itertools.product((0, 1), repeat=len(big)):
        chosen = {v for v, b in zip(big, bits) if b}
        # every edge between two unchosen big vertices is uncovered -> infeasible
        if any(u not in chosen and w in big and w not in chosen for u in big for w in adj[u]):
            continue
        # remaining vertices: those adjacent to an unchosen big vertex must be chosen
        must = {w for u in big if u not in chosen for w in adj[u]}
        rest = [v for v in adj if v not in big]
        seen = set(); total = len(chosen)
        for s in rest:
            if s in seen: continue
            comp = []; stack = [s]; seen.add(s)
            while stack:
                v = stack.pop(); comp.append(v)
                for w in adj[v]:
                    if w not in big and w not in seen:
                        seen.add(w); stack.append(w)
            total += path_vc(comp, adj, set(big), must)
        best = total if best is None else min(best, total)
    return best

def path_vc(comp, adj, big, must):
    # comp induces a path or cycle (max degree 2 within comp); exact DP over an ordering
    inner = {v: [w for w in adj[v] if w not in big] for v in comp}
    ends = [v for v in comp if len(inner[v]) <= 1]
    order = []
    start = ends[0] if ends else comp[0]
    prev = None; v = start
    while True:
        order.append(v)
        nxt = [w for w in inner[v] if w != prev and w not in order]
        if not nxt: break
        prev, v = v, nxt[0]
    assert len(order) == len(comp)
    cyc = not ends
    best = None
    for first in (0, 1):
        dp = {first: first}
        if order[0] in must and first == 0: continue
        for v in order[1:]:
            nd = {}
            for b in (0, 1):
                if v in must and b == 0: continue
                cand = [dp[a] + b for a in dp if a or b]
                if cand: nd[b] = min(cand)
            dp = nd
        for last, val in dp.items():
            if cyc and not (last or first): continue
            best = val if best is None else min(best, val)
    return best

fn = sys.argv[1]
l, d = load(fn)
assert len(d) == (l + 1) ** 2 and all(x % 2 == 0 and y % 2 == 0 for x, y in d), "every even node a parasite"
forced = set()
for p, dem in d.items():
    cs = cells_of(p, l)
    assert 1 <= dem <= len(cs)
    if dem == len(cs): forced |= set(cs)
adj = {}
for p, dem in d.items():
    cs = cells_of(p, l)
    rest = [c for c in cs if c not in forced]
    r = dem - (len(cs) - len(rest))
    if r <= 0: continue
    assert r == 1 and len(rest) == 2, ("unexpected residual constraint", p, r, rest)
    a, b = rest
    adj.setdefault(a, set()).add(b); adj.setdefault(b, set()).add(a)
k = vc_exact(adj)
print(fn.split('/')[-1], "forced", len(forced), "residual graph", len(adj), "vertices; min emitters =", len(forced) + k)

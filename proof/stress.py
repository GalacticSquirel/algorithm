"""Stress test for the gadgets: random orthogonal drawings at the README's scale.

The README says the gadget invariants (I1)-(I3) hold for every drawing, not just the
hand-made ones in gen.GRAPHS. This builds random drawings of planar graphs of maximum
degree 3 (random subgraphs of a k x k grid, degree-4 points thinned, every cycle
given a vertex, no loops; parallel edges allowed), runs gen.instance on each (which
asserts (I1)-(I3)), then rereads the instance the way check.py does and checks:
- the forced cells are exactly the non-choice cells (count F);
- every residual constraint is "one of these two cells";
- the residual graph is G with each edge replaced by a path with an even number of
  inner vertices, edge for edge.

Run: python3 stress.py [seed] [count]
Not part of the Java submission."""
import random, sys
import gen


def drawing(rnd, k):
    E = set()
    for x in range(k):
        for y in range(k):
            for d in ((1, 0), (0, 1)):
                if x + d[0] < k and y + d[1] < k and rnd.random() < 0.55:
                    E.add(((x, y), (x + d[0], y + d[1])))
    adj = {}
    for a, b in E:
        adj.setdefault(a, set()).add(b)
        adj.setdefault(b, set()).add(a)
    for v in list(adj):
        while len(adj[v]) > 3:
            w = rnd.choice(sorted(adj[v]))
            adj[v].discard(w)
            adj[w].discard(v)
    adj = {v: s for v, s in adj.items() if s}
    verts = {v for v in adj if len(adj[v]) != 2 or rnd.random() < 0.25}

    def trace(u, w):
        p = [u, w]
        while p[-1] not in verts:
            p.append(next(z for z in adj[p[-1]] if z != p[-2]))
        return p

    while True:  # give every cycle a vertex and break loops (u back to u)
        bad = None
        for v in adj:
            if v in verts:
                continue
            seen, st = {v}, [v]
            while st and bad is None:
                a = st.pop()
                if a in verts:
                    break
                for b in adj[a]:
                    if b not in seen:
                        seen.add(b)
                        st.append(b)
            else:
                bad = v
            if bad:
                break
        if bad is None:
            for u in verts:
                for w in adj[u]:
                    p = trace(u, w)
                    if p[-1] == u:
                        bad = p[len(p) // 2]
                        break
                if bad:
                    break
        if bad is None:
            break
        verts.add(bad)
    name = lambda p: f"v{p[0]}_{p[1]}"
    vs = {name(v): v for v in verts}
    es, done = [], set()
    for u in verts:
        for w in adj[u]:
            p = trace(u, w)
            key = frozenset([(p[0], p[1]), (p[-1], p[-2])])
            if key not in done:
                done.add(key)
                es.append((name(u), name(p[-1]), p))
    return vs, es


def residual(l, demand):
    def block(p):
        i, j = p
        return [(a, b) for a in (i - 1, i) for b in (j - 1, j) if 0 <= a < l and 0 <= b < l]
    forced = set()
    for p, d in demand.items():
        cs = block(p)
        assert 1 <= d <= len(cs)
        if d == len(cs):
            forced |= set(cs)
    adj = {}
    for p, d in demand.items():
        cs = block(p)
        rest = [c for c in cs if c not in forced]
        r = d - (len(cs) - len(rest))
        if r > 0:
            assert r == 1 and len(rest) == 2, (p, r, rest)
            a, b = rest
            adj.setdefault(a, set()).add(b)
            adj.setdefault(b, set()).add(a)
    return forced, adj


def check(vs, es):
    l, demand, F, s, V, _ = gen.instance(vs, es)
    forced, adj = residual(l, demand)
    assert len(forced) == F
    vcell = {v: (gen.S * x + gen.MARGIN, gen.S * y + gen.MARGIN) for v, (x, y) in vs.items()}
    c2v = {c: v for v, c in vcell.items()}
    assert all(len(adj[c]) == 2 for c in adj if c not in c2v)
    got = []
    for v, c in vcell.items():
        assert len(adj.get(c, ())) == sum((e[0] == v) + (e[1] == v) for e in es)
        for w in adj.get(c, ()):
            p = [c, w]
            while p[-1] not in c2v:
                p.append(next(x for x in adj[p[-1]] if x != p[-2]))
            assert (len(p) - 2) % 2 == 0
            got.append(tuple(sorted((v, c2v[p[-1]]))))
    assert sorted(got) == sorted(tuple(sorted(e[:2])) for e in es for _ in (0, 1))


if __name__ == "__main__":
    seed = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    count = int(sys.argv[2]) if len(sys.argv) > 2 else 300
    rnd = random.Random(seed)
    n = 0
    for _ in range(count):
        vs, es = drawing(rnd, rnd.randint(3, 7))
        if es:
            check(vs, es)
            n += 1
    print(f"scale {gen.S}: {n} random drawings, all gadget invariants and residual graphs correct")

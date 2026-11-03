# CS31920 assignment: emitter placement

Java solution for the CS31920 Advanced Algorithms assignment (2026-27). It reads
energy grids with parasites and prints the fewest emitters that neutralise every
parasite, or `no solution`.

## Build and run

Needs Java 21 or later; the project targets Java 26.

```
javac -d out src/uk/ac/aber/cs31920/assignment/*.java
java -cp out uk.ac.aber.cs31920.assignment.MySolution inputs/input01.txt inputs/input02.txt
```

Each input file gives one output line. Add `-Dverbose=1` before `-cp` to print the
runtime, emitter count, lower bound and whether the answer is proven optimal to
standard error.

Only the Java Collections Framework is used.

A file that cannot be opened, or is not in the input format (first line `l`, then
`x y demand` lines), prints `no solution`. A grid with no parasites prints an empty
line (zero emitters). If a node is listed twice, the later line counts.

## How it works

| Class | Role |
|---|---|
| `MySolution` | Entry point: one line of output per input file |
| `Instance` | Reads an input file |
| `Grid` | Triangle mesh: node numbering, coordinates, neighbours |
| `FlowNetwork` | Dinitz max flow, plus Ford-Fulkerson searches from one changed arc and an undo log |
| `PulseNetwork` | Node-split flow network that checks a set of emitters |
| `EmitterPlacer` | No-solution test, forced emitters, lower bounds, greedy, local search, branch and bound |
| `RelocationSearch` | Iterated local search for a placement with one emitter fewer |
| `RegionBound` | Lower bound from small regions that must hold an emitter |

1. **No-solution test (exact).** The instance is solvable exactly when every
   parasite's demand is at most its number of parasite-free neighbours.
2. **Checking a placement.** Pulses are vertex-disjoint paths, so a placement works
   exactly when the maximum flow in the node-split network equals the total demand
   (22 Oct slides 189-191, 195). From scratch this uses Dinitz's algorithm (28 Oct
   slides 229-234). After adding one emitter, only augmenting paths that start at it
   can exist (cut argument, 23 Oct slides 215-219), so Ford-Fulkerson (22 Oct slide
   203) searches from that emitter only. After removing one, paths can only end at
   the parasites that lost a pulse, so the search runs backwards from them.
3. **Forced emitters and lower bounds.** A free node next to two parasites whose
   demand equals their free degree must be an emitter. After placing them, each
   parasite that is still short has a set of "helpers": the free nodes that can
   reach it in the residual network without passing through the sink. By the
   max-flow min-cut argument (23 Oct slides 215-219) every solution adds an emitter
   in each helper set. So the lower bound is the forced emitters plus the larger of
   ceil(missing pulses / 8) and the number of pairwise disjoint helper sets. A
   stronger bound is used when the
   greedy answer does not meet it (a relaxation, 16 Oct slides 148 and 153, used as
   the lower bound module of branch and bound, 14 Oct slide 117):
   - region bound: a relaxed max flow, where every free node outside a small region
     may send pulses, shows when that region must hold an emitter. Near crowded
     corners and edges such a region has only low-degree nodes, so its emitter
     wastes capacity.
4. **Upper bound.** Greedy construction plus a local search that drops unneeded
   emitters (14 Oct slide 110). The greedy ranks candidates with a cheap local
   estimate (tried with an undo log, then rolled back) and checks only the best few
   exactly. An iterated local search (14 Oct slides 109-110) then moves emitters
   around to find a placement with one emitter fewer, as long as it can.
5. **Branch and bound** (14 Oct slides 117-119), seeded with the greedy result. At
   each node it branches on the helpers of the short parasite with the fewest
   helpers, since every completion contains one of them (same min-cut argument),
   which loses no solution. A branch is cut off when the emitters left are fewer
   than the bound from step 3, or when the nodes ruled out so far leave too little
   room even in the relaxed network.
   It always runs until the answer
   is proven minimal, as the brief requires. For experiments, `-Dbudget=1000` caps
   local search and branch and bound at 1000 ms per input; the best valid placement
   found so far is then printed, without a guarantee that it is minimal.

## Complexity

N is the number of grid nodes ((l+1)^2 + l^2), so the flow network has O(N) nodes
and O(N) arcs. D is the total demand (at most 8 per parasite), P the number of
parasites, K the number of emitters and M = 40 the greedy candidates per step.

| Step | Before | After |
|---|---|---|
| Max flow from scratch | Ford-Fulkerson O(D·N) (23 Oct s212), O(N^2) when D grows with N | Dinitz: O(n^2 e) in general (28 Oct s234), O(sqrt(N)·N) on this near-unit network |
| Add one emitter | copy O(N) + search from the source O(N) | search from the new emitter only: O(N) worst case, usually just the area around it |
| Greedy step | M exact trials, O(M·N) | one O(N) reachability search + M estimates over the reachable area R, O(N + M·R), with undo instead of copying + 4 exact trials |
| Local search pass | rebuild after every success, each check a full max flow: O(K·D·N) per check round | one network kept up to date; each check is a backward search from at most 8 parasites, O(N), skipped outright when the others lack spare capacity |
| Branch and bound node overhead | O(N·P) parasite distances + sort | O(N log N) sort (distances precomputed once in O(N)) |
| Whole minimisation | exponential worst case | exponential worst case; no polynomial exact method exists unless P = NP (see below) |

The sqrt(N) phase bound for Dinitz is not on the slides. The argument is in the
`FlowNetwork` class comment: free nodes carry at most one unit and other nodes at
most 8, so after d phases at most 8N/d units are missing.

Measured on generated grids with sparse parasites (total time per file):

| Grid | Before | After |
|---|---|---|
| l = 40 (150 parasites) | 0.3-0.6 s | 0.15 s |
| l = 80 (600 parasites) | 5.8-7.1 s | 0.8 s |
| l = 120 (1500 parasites) | 43-45 s | 4.3 s |

On the ten given inputs both versions take 30-90 ms of solver time per file, which is
mostly JVM start-up.

## Results on the given inputs

| Input | Output | Emitters | Lower bound | Proven optimal |
|---|---|---|---|---|
| input01 | `4 10 2 2` | 2 | 2 | yes |
| input02 | `0 2 4 2` | 2 | 2 | yes |
| input03 | `18 18` | 1 | 1 | yes |
| input04 | 34 emitters | 34 | 34 | yes |
| input05 | `no solution` | | | |
| input06 | `10 6 14 8 18 14 4 2` | 4 | 4 | yes |
| input07 | 9 emitters | 9 | 9 | yes |
| input08 | `no solution` | | | |
| input09 | 9 emitters | 9 | 9 | yes |
| input10 | 7 emitters | 7 | 7 | yes |

## Why these answers are minimal

Each answer below is matched by a lower bound, so no placement with fewer emitters
exists. One emitter on node v sends at most degree(v) pulses, so with total demand D
and largest free degree d at least ceil(D / d) emitters are needed.

| Input | Proof |
|---|---|
| input01 | D = 11, d = 8, so at least 2 |
| input02 | D = 7. The only degree-8 node, (2,2), holds a parasite, so d = 5 and at least 2 |
| input03 | D = 6 > 0, so at least 1 |
| input04 | 32 free nodes each touch two parasites whose demand equals their free degree, so all 32 are emitters in every solution. Those 32 leave 3 pulses missing. Parasite (2,14) can then only be helped from the corner pocket (0,14), (0,16), (1,15), (2,16), and parasite (10,6) only from nodes outside that pocket. The two helper sets are disjoint, so at least 32 + 2 = 34 |
| input05 | Parasite (10,10) needs 5 pulses but has only 4 neighbours without a parasite. Every pulse enters through a different one of them, so no placement works |
| input06 | D = 31, d = 8, so at least 4 |
| input07 | D = 68, d = 8, so at least 9 |
| input08 | Parasite (24,26) needs 6 pulses but has only 5 neighbours without a parasite, so no placement works |
| input09 | D = 71, d = 8, so at least 9 |
| input10 | D = 56, d = 8, so at least 7 |

The program checks every placement it prints with an exact maximum flow, so each
answer works and is the minimum.

## Finding the minimum is NP-hard

**Theorem.** The decision problem EMITTERS asks: given an instance and a number K,
can at most K emitters neutralise every parasite? EMITTERS is NP-hard, and it is
NP-complete when the size of an instance is measured by its grid size
N = (l+1)^2 + l^2. So no algorithm finds the smallest number of emitters in
polynomial time unless P = NP.

*Which input size.* The two halves of the theorem need different care.
- Hardness holds under any reasonable encoding. The instances built below list every
  even node as a parasite, so their files have length at least l^2, and every number
  in them (l, the coordinates, demands of at most 4) is polynomial in the size of the
  vertex cover instance. So the reduction is polynomial whether l is written in
  binary or in unary. (This also makes the problem *strongly* NP-hard: it stays hard
  when all numbers are small.)
- Membership in NP is proved only for the grid-size measure, which is the same as
  writing l in unary. If l is written in binary, a file of a few lines can describe a
  grid exponentially larger than the file, and the verifier below, which builds the
  whole grid, is then not polynomial. We do not claim membership for that encoding,
  and the hardness result does not need it.

Definitions used: NP by certificates (2 Oct, slide 33), NP-hardness (2 Oct, slide
38), and proving NP-hardness by a polynomial reduction from a known NP-hard problem,
since reductions are transitive (2 Oct, slides 39-40). The problem reduced from, and
the drawing theorem, come from the literature, not the slides.

### EMITTERS is in NP

The certificate is a set of at most K free nodes (it never needs more than N). A
verifier builds the node-split network for it, which has O(N) nodes and edges, and
checks that its maximum flow equals the total demand. This is the
check `PulseNetwork` does (22 Oct, slides 189-195), and it runs in polynomial time,
for example with Dinitz's algorithm in O(n^2 e) (28 Oct, slide 234).

### Reduction from planar vertex cover

VERTEX COVER is NP-complete even for planar graphs of maximum degree 3 (Garey and
Johnson 1977; Mohar 2001 gives a complete proof even for planar graphs where every
vertex has degree exactly 3). Degree at most 3 is needed below: a vertex cell with
four arms would leave no designated corner for its arms. Let G be such a graph and k a number. The construction below runs in
polynomial time.

**Cells and corners.** Call the odd nodes *cells*: (2i+1, 2j+1) is the centre of a
unit square. Call the even nodes *corners*. A cell is adjacent to exactly its four
corners and to nothing else. A corner is adjacent to the (up to) four cells around it
and to other corners.

**Lemma 1.** Suppose every even node holds a parasite and every odd node is free.
Then a set S of cells works if and only if every corner p has at least d(p) cells of
S among the cells around it.

*Proof.* Emitters can only stand on free nodes, which are cells. A pulse leaves an
emitter cell along an edge into a corner. That corner holds a parasite, so a pulse
cannot pass through it, and the pulse ends there. So every pulse is a single edge,
and corner p receives at most one pulse from each cell around it, since there is at
most one pulse per edge. Conversely, if p has at least d(p) emitter cells around it,
let each of them send one pulse along its edge to p. Distinct corners use distinct
edges of an emitter, and single-edge pulses share no node except their ends, so all
the rules hold. ∎

**Lemma 2.** In such an instance, if d(p) equals the number of cells around p, then
every cell around p is an emitter in every solution. (This follows directly from
Lemma 1.)

**Lemma 3.** Replace an edge uv of a graph H by a path u-a-b-v. The minimum vertex
cover grows by exactly 1.

*Proof.* Take a cover C of H. If u is in C, then C + b covers the new graph. If v is
in C, then C + a does.

Conversely, take a cover C' of the new graph. It contains a or b, because of the edge
ab.
- If C' contains u or v, then C' without a and b covers H, and it is at least 1
  smaller.
- If it contains neither u nor v, it must contain both a and b, because of the edges
  ua and bv. Then C' − {a, b} + {u} covers H, and it is exactly 1 smaller. ∎

Applying Lemma 3 t times to one edge replaces that edge by a path with 2t inner
vertices, and the minimum cover grows by t.

**Construction.**

1. Draw G orthogonally on a grid. Vertices go to grid points and edges to
   non-crossing rectilinear paths. This takes polynomial time and O(n^2) area for
   planar graphs of maximum degree 4 (Valiant 1981; Tamassia 1987). Scale the drawing
   by 32, so every straight segment is at least 32 cells long and separate parts of
   the drawing are at least 32 cells apart. The gadgets below need at most 4 cells
   next to a vertex and a straight stretch of 20 cells, and reach at most 3 cells
   sideways, so they always fit and only touch other chains where chains meet at a
   vertex. (`gen.py` raises an error if a gadget does not fit. The
   cube graph Q3 needs a scale of at least 27.)
2. Mark *choice cells*. Each vertex becomes one cell, and each edge becomes a chain of
   cells along its path, consecutive cells being neighbours. Two local changes are
   needed:
   - **T-junction.** A degree-3 vertex has two collinear arms and one perpendicular
     arm. One collinear arm steps one row away from the perpendicular arm by a
     diagonal step, runs two cells, and steps back.
   - **Parity fix.** If an edge's chain has an odd number of inner cells, one straight
     run takes a detour: three diagonal steps out, two straight cells, then three
     cells straight back. That adds 3 cells, so the inner count becomes even.

   ```
   T-junction (V = vertex cell)        parity fix on a run going right
   . . # . . . . .                     . . . # # # # .
   . . # . . . . .                     . . # . . . # .
   # # V . . . # #                     . # . . . . # .
   . . . # # # . .                     # . . . . . # #
   ```
3. **Demands.** For every corner p, let f(p) be the number of non-choice cells around
   it. For each pair of consecutive cells a, b in a chain, choose one *designated*
   corner whose choice cells are exactly {a, b}. Set d(p) = f(p) + 1 at designated
   corners and d(p) = f(p) everywhere else. Every even node becomes a parasite with
   demand d(p), and every odd node stays free.

The gadgets keep three invariants. Each can be checked locally on the pictures above
and on straight runs and turns. `proof/gen.py` asserts all three on every instance it
builds, and `proof/stress.py` builds 1,200 random drawings to check that they hold
beyond the hand-made examples.

- **(I1)** No 2×2 block consists of choice cells only. So f(p) ≥ 1 and 1 ≤ d(p) ≤ 4
  at every corner, which makes the instance valid, and it passes the no-solution test.
- **(I2)** Every non-choice cell lies in a 2×2 block, or a board-edge block, of
  non-choice cells. That block's corner has d = f = its number of cells, so the cell
  is forced (Lemma 2).
- **(I3)** Every pair of consecutive chain cells has a designated corner.

**Correctness.** By (I2) and Lemma 2, all non-choice cells, say F of them, are
emitters in every solution. With them in place:
- a corner with d(p) = f(p) is already satisfied;
- a designated corner still needs one more emitter, namely a or b.

So, by Lemma 1, a set of emitters works exactly when it contains F and a vertex cover
of G'. Here G' is G with every edge e replaced by a path with 2t_e inner vertices,
and t_e is half the number of inner cells of e's chain. By Lemma 3:

    minimum number of emitters = F + VC(G) + Σ t_e.

So G has a vertex cover of size k if and only if the instance can be solved with
F + k + Σ t_e emitters. The instance has O(n^2) nodes and is built in polynomial
time. Together with membership in NP, this proves the theorem. ∎

**Checked by program.** The proof above is the argument. These scripts look for
mistakes in it; they are not part of the Java submission.
- `proof/gen.py` builds instances for seven graphs: K4, a 5-cycle, a small tree, the cube Q3 (eight T-junctions, with the perpendicular arm in
  all four directions), K2,3, the triangular prism, and a path with several bends next
  to an isolated vertex. It asserts (I1)-(I3) on each.
- `proof/check.py` reads only the instance file, rederives the forced cells and the
  residual pair constraints, and computes the minimum exactly. It matches
  F + VC(G) + Σ t_e for all seven.
- `proof/stress.py` does the same on 1,200 random orthogonal drawings (seeds 0-3, 300
  each), including parallel edges, bends next to vertices and isolated vertices. It
  also checks that the residual graph is G with each edge replaced by an even-length
  path, edge for edge, not just that the numbers agree. The test was checked to fail
  when the T-junction step or the designated corners are removed from `gen.py`.
- `check.py` and `stress.py` rely on Lemmas 1 and 2. `proof/semantic.py` does not. It
  has its own node-split max flow written from the brief's rules (vertex-disjoint
  pulses, one per incident edge, nothing through an emitter or a parasite). First it
  checks Lemma 1 against max flow on 300 random small grids. Then, at scale 14, for
  six of the graphs, it finds by max flow which cells are needed in every solution
  and which pairs of cells must contain an emitter, checks that the pairs form G with
  even subdivisions, and proves the minimum from the rules alone: a lower bound from
  those facts plus an exact vertex cover, and an emitter set of that size that passes
  max flow. Each minimum equals F + VC(G) + Σ t_e.
- The Java solver, an independent max-flow implementation, proves the same minimum on
  three scale-5 instances (one edge 193, a path of three vertices 190, a 4-cycle 186;
  build them with `gen.S = 5`).
- On the full-size K4 instance the solver did not finish within ten minutes. That
  fits the theorem but is not evidence for it.

**References.**
- M. R. Garey and D. S. Johnson, "The rectilinear Steiner tree problem is
  NP-complete", SIAM J. Appl. Math. 32(4), 826-834, 1977 (vertex cover in planar
  graphs of degree at most 3). Do not cite Garey, Johnson and Stockmeyer 1976 for
  this: their planar vertex cover reduction only reaches maximum degree 6.
- B. Mohar, "Face covers and the genus problem for apex graphs", J. Combin. Theory
  Ser. B 82(1), 102-117, 2001 (a complete proof for planar graphs with every degree
  exactly 3, after Theorem 4.1).
- L. G. Valiant, "Universality considerations in VLSI circuits", IEEE Trans.
  Computers C-30(2), 1981.
- R. Tamassia, "On embedding a graph in the grid with the minimum number of bends",
  SIAM J. Comput. 16(3), 1987.

## AI use

This code was written with Claude (Anthropic), an AI assistant, at Owen's request.
This has to be disclosed in the report's References section, with the significant
prompts and how much of the output was used.

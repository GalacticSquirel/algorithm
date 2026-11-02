# CS31920 assignment: emitter placement

Java solution for the CS31920 Advanced Algorithms assignment (2026-27). It reads
energy grids with parasites and prints the fewest emitters that neutralise every
parasite, or `no solution`.

## Build and run

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
| Whole minimisation | exponential worst case | exponential worst case (no polynomial exact method is known) |

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

## AI use

This code was written with Claude (Anthropic), an AI assistant, at Owen's request.
This has to be disclosed in the report's References section, with the significant
prompts and how much of the output was used.

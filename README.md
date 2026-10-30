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

## How it works

| Class | Role |
|---|---|
| `MySolution` | Entry point: one line of output per input file |
| `Instance` | Reads an input file |
| `Grid` | Triangle mesh: node numbering, coordinates, neighbours |
| `FlowNetwork` | Ford-Fulkerson max flow with BFS augmenting paths |
| `PulseNetwork` | Node-split flow network that checks a set of emitters |
| `EmitterPlacer` | No-solution test, forced emitters, lower bound, greedy, local search, branch and bound |

1. **No-solution test (exact).** The instance is solvable exactly when every
   parasite's demand is at most its number of parasite-free neighbours.
2. **Checking a placement.** Pulses are vertex-disjoint paths, so a placement works
   exactly when the maximum flow in the node-split network equals the total demand
   (22 Oct slides 189-191, 195, 203; 23 Oct slide 212).
3. **Forced emitters and lower bound.** A free node next to two parasites whose
   demand equals their free degree must be an emitter. The lower bound is the
   forced emitters plus ceil(missing pulses / 8).
4. **Upper bound.** Greedy construction plus a local search that drops unneeded
   emitters (14 Oct slide 110).
5. **Branch and bound** (14 Oct slides 117-119), seeded with the greedy result. It
   branches only on nodes on the sink side of the minimum cut (max-flow min-cut,
   23 Oct slides 215-219), which loses no solution. It has a 1 second budget per
   input; if that runs out the best valid placement found so far is printed.

## Results on the given inputs

| Input | Output | Emitters | Lower bound | Proven optimal |
|---|---|---|---|---|
| input01 | `4 10 2 2` | 2 | 2 | yes |
| input02 | `0 2 4 2` | 2 | 2 | yes |
| input03 | `18 18` | 1 | 1 | yes |
| input04 | 34 emitters | 34 | 33 | yes (search for 33 finished) |
| input05 | `no solution` | | | |
| input06 | `10 6 14 8 18 14 4 2` | 4 | 4 | yes |
| input07 | 9 emitters | 9 | 9 | yes |
| input08 | `no solution` | | | |
| input09 | 9 emitters | 9 | 9 | yes |
| input10 | 7 emitters | 7 | 7 | yes |

## AI use

This code was written with Claude (Anthropic), an AI assistant, at Owen's request.
This has to be disclosed in the report's References section, with the significant
prompts and how much of the output was used.

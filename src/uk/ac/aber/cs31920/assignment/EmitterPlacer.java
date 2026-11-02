package uk.ac.aber.cs31920.assignment;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

/**
 * Finds a placement with as few emitters as possible.
 *
 * The steps are:
 *
 * 1. Exact no-solution test. If every free node held an emitter, each parasite
 *    could take one pulse straight from every neighbour without a parasite,
 *    and these one-edge pulses never share a node. So the instance is
 *    solvable if and only if every parasite's demand is at most its number of
 *    parasite-free neighbours. (If a parasite fails this, no placement can
 *    work, because each pulse enters it through a different edge from a
 *    non-parasite node.)
 * 2. Forced emitters and a lower bound: some free nodes must hold an
 *    emitter in every solution (see forcedEmitters). One emitter sends at most
 *    degree(v) <= 8 pulses, which gives the lower bound (see lowerBound).
 * 3. Upper bound: starting from the forced emitters, a greedy construction
 *    keeps adding the emitter that delivers the most extra pulses (ranked
 *    with a quick local estimate, then checked exactly). A local search
 *    (14 Oct, slide 110) then removes emitters that are no longer needed,
 *    and an iterated local search (14 Oct, slides 109-110, see
 *    RelocationSearch) moves emitters around to look for a placement with
 *    one emitter fewer.
 * 4. If the best placement is still above the lower bound, the bound is
 *    raised with a relaxation (see RegionBound).
 * 5. Branch and bound (14 Oct, slides 117-118), seeded with the best
 *    placement so far as the incumbent (slide 119): it keeps asking for a
 *    working set with one emitter fewer. A branch is cut off when the pulses
 *    still missing need more emitters than are left, or when the nodes ruled
 *    out so far leave too little room even in a relaxed network. Branching
 *    happens only on nodes allowed by the minimum cut argument in
 *    PulseNetwork.usefulNextEmitters, so the cut-down choice loses no
 *    solution. When a search ends without finding a smaller set, or the best
 *    set meets the lower bound, the answer is proven optimal.
 *
 * Every candidate set is checked with the maximum flow in PulseNetwork, so any
 * placement this class returns is guaranteed to work.
 *
 * The brief asks for the smallest possible number of emitters, so by default
 * the search runs until the answer is proven optimal, however long that takes
 * (it is exponential in the worst case). An optional time budget can be set
 * for experiments; if it runs out, the best placement found so far is
 * returned, which is valid but not proven optimal.
 */
public final class EmitterPlacer {

    /** Candidates estimated per greedy step. */
    private static final int GREEDY_CANDIDATES = 40;
    /** Of those, how many of the best estimates are then checked exactly. */
    private static final int EXACT_CHECKS = 4;
    /** Moves the iterated local search may try for each smaller size. */
    private static final int RELOCATION_MOVES = 3000;

    private final Instance instance;
    private final Grid grid;
    private final long budgetMillis;
    /** The most pulses one emitter can send: the largest degree of a free node (8 on most grids). */
    private final int maxEmitterPulses;
    /** parasiteDistance[v]: distance from v to the nearest parasite (for ordering only). */
    private final int[] parasiteDistance;
    // Set by solve().
    private long deadline;
    private boolean timedOut;
    private boolean provenOptimal;
    private List<Integer> forced = new ArrayList<>();
    private int lower;

    /**
     * @param instance     the instance to solve
     * @param budgetMillis time allowed for local search and branch and bound,
     *                     or Long.MAX_VALUE for no limit
     */
    public EmitterPlacer(Instance instance, long budgetMillis) {
        this.instance = instance;
        this.grid = instance.grid();
        this.budgetMillis = budgetMillis;
        int best = 0;
        for (int v = 0; v < grid.nodeCount(); v++) {
            if (!instance.hasParasite(v)) {
                best = Math.max(best, grid.degree(v));
            }
        }
        // At least 1, so the divisions below are safe on a grid with no edges.
        this.maxEmitterPulses = Math.max(1, best);
        this.parasiteDistance = distancesToParasites();
    }

    /**
     * Solves the instance.
     *
     * @return the emitter nodes, or null if the instance has no solution
     */
    public List<Integer> solve() {
        timedOut = false;
        if (!isSolvable()) {
            return null;
        }
        forced = forcedEmitters();
        PulseNetwork start = build(forced);
        lower = Math.max(emittersNeeded(instance.totalDemand()),
            forced.size() + emittersNeeded(start.deficit()));

        List<Integer> greedy = greedy(start);
        // The time budget covers the improvement steps that follow.
        long now = System.currentTimeMillis();
        deadline = budgetMillis >= Long.MAX_VALUE - now ? Long.MAX_VALUE : now + budgetMillis;
        List<Integer> best = localSearch(greedy);

        // Iterated local search: while it keeps finding a working placement
        // with one emitter fewer, take it.
        RelocationSearch relocation = new RelocationSearch(instance, forced);
        while (best.size() > lower && System.currentTimeMillis() <= deadline) {
            List<Integer> smaller = relocation.search(best, RELOCATION_MOVES);
            if (smaller == null) {
                break;
            }
            best = smaller;
        }
        if (best.size() > lower) {
            lower = Math.max(lower, new RegionBound(instance, maxEmitterPulses).lowerBound(forced));
        }

        // Branch and bound: keep asking for a placement with one emitter fewer
        // than the best so far. When a search finishes without finding one,
        // the best placement is optimal.
        provenOptimal = best.size() == lower;
        while (!provenOptimal && !timedOut) {
            List<Integer> smaller = branchAndBound(start, best.size() - 1);
            if (smaller != null) {
                best = smaller;
                provenOptimal = best.size() == lower;
            } else if (!timedOut) {
                provenOptimal = true;
            }
        }
        return best;
    }

    /** True if the last call to solve proved that its answer is optimal. */
    public boolean provenOptimal() {
        return provenOptimal;
    }

    /** Step 1: the exact test described in the class comment. */
    public boolean isSolvable() {
        for (int p : instance.parasites()) {
            int freeNeighbours = 0;
            for (int w : grid.neighbours(p)) {
                if (!instance.hasParasite(w)) {
                    freeNeighbours++;
                }
            }
            if (instance.demand(p) > freeNeighbours) {
                return false;
            }
        }
        return true;
    }

    /**
     * The lower bound worked out by solve(). Step 2 gives the larger of
     * ceil(total demand / most pulses one emitter can send) and the number of
     * forced emitters plus ceil(pulses they cannot deliver / most pulses one
     * emitter can send). The second part is valid because adding an emitter v
     * raises the maximum flow by at most degree(v): removing v's own pulses
     * from any flow leaves a flow for the other emitters. Step 4 may raise it
     * further (RegionBound).
     */
    public int lowerBound() {
        return lower;
    }

    /**
     * Nodes that hold an emitter in every solution.
     *
     * Call a parasite tight if its demand equals its number of parasite-free
     * neighbours. Every such neighbour must then deliver a pulse to it, either
     * as an emitter or as the last node of a pulse. A free node next to two
     * tight parasites cannot be the last node of two pulses (it can carry only
     * one), so it has to be an emitter.
     */
    private List<Integer> forcedEmitters() {
        int n = grid.nodeCount();
        int[] tightNeighbours = new int[n];
        for (int p : instance.parasites()) {
            int free = 0;
            for (int w : grid.neighbours(p)) {
                if (!instance.hasParasite(w)) {
                    free++;
                }
            }
            if (free == instance.demand(p)) {
                for (int w : grid.neighbours(p)) {
                    tightNeighbours[w]++;
                }
            }
        }
        List<Integer> result = new ArrayList<>();
        for (int v = 0; v < n; v++) {
            if (!instance.hasParasite(v) && tightNeighbours[v] >= 2) {
                result.add(v);
            }
        }
        return result;
    }

    private int emittersNeeded(int pulses) {
        return (pulses + maxEmitterPulses - 1) / maxEmitterPulses;
    }

    /**
     * Step 3a: greedy construction. Each round looks at the free nodes allowed
     * by the minimum cut argument, ranks them by how many of their edges could
     * still carry a pulse to a parasite, tries the best few for real and keeps
     * the one that delivers the most extra pulses.
     */
    private List<Integer> greedy(PulseNetwork start) {
        PulseNetwork network = start.copy();
        List<Integer> emitters = new ArrayList<>(forced);
        while (!network.isFeasible()) {
            boolean[] reach = network.residualReachToSink();
            List<int[]> estimates = new ArrayList<>();
            int chosen = -1;
            for (int v : rankedCandidates(network, reach, null, GREEDY_CANDIDATES)) {
                // A quick estimate never overstates the gain, so an estimate
                // that uses every edge of v is already the best possible.
                int estimate = network.estimateGain(v, reach);
                if (estimate == grid.degree(v)) {
                    chosen = v;
                    break;
                }
                estimates.add(new int[] {v, estimate});
            }
            if (chosen < 0) {
                // Check the few best estimates exactly and keep the winner.
                estimates.sort(Comparator.comparingInt((int[] e) -> -e[1]));
                int bestGain = -1;
                for (int i = 0; i < Math.min(EXACT_CHECKS, estimates.size()); i++) {
                    int v = estimates.get(i)[0];
                    int gain = network.exactGain(v);
                    if (gain > bestGain) {
                        bestGain = gain;
                        chosen = v;
                    }
                }
            }
            if (chosen < 0) {
                // Cannot happen for a solvable instance: by the minimum cut
                // argument some free node can always still raise the flow.
                throw new IllegalStateException("no candidate emitter left");
            }
            // The chosen emitter is added exactly, so the flow stays maximum.
            network.addEmitter(chosen);
            emitters.add(chosen);
        }
        return emitters;
    }

    /**
     * Step 3b: local search. Tries to drop each emitter in turn and keeps the
     * smaller set if the rest still neutralise every parasite. Repeats until
     * no emitter can be dropped.
     */
    private List<Integer> localSearch(List<Integer> emitters) {
        List<Integer> current = new ArrayList<>(emitters);
        if (current.size() == lower) {
            return current; // already optimal, nothing to drop
        }
        // Build the network once and keep it in step with every change, so a
        // pass costs one removal check per emitter.
        PulseNetwork network = build(current);
        boolean improved = true;
        while (improved) {
            improved = false;
            for (int i = 0; i < current.size(); i++) {
                if (System.currentTimeMillis() > deadline) {
                    return current; // out of time; current is still valid
                }
                // Quick necessary condition: the others must have room for
                // the pulses this emitter sends now.
                int v = current.get(i);
                if (network.spareCapacityWithout(v) < network.pulsesFrom(v)) {
                    continue;
                }
                PulseNetwork trial = network.copy();
                trial.removeEmitter(current.get(i));
                if (trial.isFeasible()) {
                    network = trial;
                    current.remove(i);
                    i--;
                    improved = true;
                }
            }
        }
        return current;
    }

    /** Builds the network with the given emitters placed. */
    private PulseNetwork build(List<Integer> emitters) {
        PulseNetwork network = new PulseNetwork(instance);
        for (int v : emitters) {
            network.addEmitter(v);
        }
        return network;
    }

    /**
     * Step 5: depth-first branch and bound for a working set of at most k
     * emitters that includes the forced ones.
     *
     * @return the emitters, or null if none exists (or the time ran out)
     */
    private List<Integer> branchAndBound(PulseNetwork start, int k) {
        List<Integer> chosen = new ArrayList<>(forced);
        boolean[] forbidden = new boolean[grid.nodeCount()];
        if (search(start.copy(), k, chosen, forbidden)) {
            return chosen;
        }
        return null;
    }

    /**
     * One node of the search tree.
     *
     * @param network   the flow network with the emitters in chosen placed
     * @param k         the number of emitters allowed in total
     * @param chosen    the emitters placed so far (extended on success)
     * @param forbidden nodes that earlier siblings already covered; a node
     *                  explored in one branch is excluded from the branches
     *                  after it, so every set is looked at only once
     */
    private boolean search(PulseNetwork network, int k, List<Integer> chosen, boolean[] forbidden) {
        if (network.isFeasible()) {
            return true;
        }
        if (System.currentTimeMillis() > deadline) {
            timedOut = true;
            return false;
        }
        int left = k - chosen.size();
        int deficit = network.deficit();
        if (emittersNeeded(deficit) > left) {
            return false;
        }

        // Exact gain of every allowed candidate, best first.
        List<int[]> children = new ArrayList<>();
        for (int v : rankedCandidates(network, network.residualReachToSink(), forbidden, Integer.MAX_VALUE)) {
            children.add(new int[] {v, network.exactGain(v)});
        }
        children.sort(Comparator.comparingInt((int[] c) -> -c[1]));

        List<Integer> excluded = new ArrayList<>();
        boolean found = false;
        for (int i = 0; i < children.size(); i++) {
            // Bound: after child i, each of the other left - 1 emitters adds
            // at most maxEmitterPulses (removing an emitter's own pulses from
            // a flow leaves a valid flow for the rest). The gains are sorted,
            // so once a child falls short every later child does too.
            if (children.get(i)[1] + (left - 1) * maxEmitterPulses < deficit) {
                break;
            }
            // Every node forbidden so far stays empty in this branch. If the
            // parasites cannot be served even with every other free node
            // allowed to send pulses, neither this branch nor any later one
            // (which forbid even more) can work.
            if (i > 0 && !coverableWithout(forbidden)) {
                break;
            }
            int v = children.get(i)[0];
            PulseNetwork next = network.copy();
            next.addEmitter(v);
            chosen.add(v);
            if (search(next, k, chosen, forbidden)) {
                found = true;
                break;
            }
            chosen.remove(chosen.size() - 1);
            if (timedOut) {
                break;
            }
            forbidden[v] = true;
            excluded.add(v);
        }
        // Undo the exclusions so that other parts of the tree are unaffected.
        for (int v : excluded) {
            forbidden[v] = false;
        }
        return found;
    }

    /**
     * Relaxation check (16 Oct, slides 148 and 153): can all demand be met if
     * every free node that is not forbidden may send pulses (and still let
     * others pass)? Any placement that avoids the forbidden nodes gives a
     * flow in this relaxed network, so if its maximum flow is below the total
     * demand, no such placement exists. Costs one maximum flow (Dinitz,
     * 28 Oct, slides 229-234).
     */
    private boolean coverableWithout(boolean[] forbidden) {
        PulseNetwork relaxed = new PulseNetwork(instance);
        for (int v = 0; v < grid.nodeCount(); v++) {
            if (!instance.hasParasite(v) && !forbidden[v]) {
                relaxed.openRelaxedSource(v);
            }
        }
        return relaxed.maximiseRelaxed() >= instance.totalDemand();
    }

    /**
     * Free nodes allowed by the minimum cut argument, ordered so the most
     * promising come first: more edges that could carry a new pulse, then
     * closer to a parasite. Returns at most limit nodes.
     */
    private List<Integer> rankedCandidates(PulseNetwork network, boolean[] reach, boolean[] forbidden,
            int limit) {
        boolean[] useful = network.usefulNextEmitters(reach);
        List<Integer> candidates = new ArrayList<>();
        int n = grid.nodeCount();
        int[] score = new int[n];
        for (int v = 0; v < n; v++) {
            if (useful[v] && (forbidden == null || !forbidden[v])) {
                candidates.add(v);
                score[v] = network.promisingEdges(v, reach);
            }
        }
        Collections.sort(candidates, Comparator
            .comparingInt((Integer v) -> -score[v])
            .thenComparingInt(v -> parasiteDistance[v])
            .thenComparingInt(v -> v));
        if (candidates.size() > limit) {
            return new ArrayList<>(candidates.subList(0, limit));
        }
        return candidates;
    }

    /**
     * Distance max(|dx|, |dy|) from every node to its nearest parasite, in
     * coordinate units. Computed once by a breadth-first search from all
     * parasites at once over the (2l+1) x (2l+1) coordinate square with king
     * moves, which takes O(N) time, instead of comparing every candidate with
     * every parasite in every search node (O(N * P) each time).
     */
    private int[] distancesToParasites() {
        int size = 2 * grid.sideLength() + 1;
        int[] square = new int[size * size];
        Arrays.fill(square, Integer.MAX_VALUE);
        ArrayDeque<Integer> queue = new ArrayDeque<>();
        for (int p : instance.parasites()) {
            int cell = grid.x(p) * size + grid.y(p);
            square[cell] = 0;
            queue.add(cell);
        }
        while (!queue.isEmpty()) {
            int cell = queue.poll();
            int x = cell / size;
            int y = cell % size;
            for (int dx = -1; dx <= 1; dx++) {
                for (int dy = -1; dy <= 1; dy++) {
                    int nx = x + dx;
                    int ny = y + dy;
                    if (nx >= 0 && ny >= 0 && nx < size && ny < size
                            && square[nx * size + ny] == Integer.MAX_VALUE) {
                        square[nx * size + ny] = square[cell] + 1;
                        queue.add(nx * size + ny);
                    }
                }
            }
        }
        int[] distance = new int[grid.nodeCount()];
        for (int v = 0; v < distance.length; v++) {
            distance[v] = square[grid.x(v) * size + grid.y(v)];
        }
        return distance;
    }
}

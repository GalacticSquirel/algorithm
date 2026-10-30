package uk.ac.aber.cs31920.assignment;

import java.util.ArrayList;
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
 *    keeps adding the emitter that delivers the most extra pulses, followed
 *    by a local search (14 Oct, slide 110) that removes emitters that are no
 *    longer needed.
 * 4. Branch and bound (14 Oct, slides 117-118), seeded with the greedy result
 *    as the incumbent (slide 119): it keeps asking for a working set with one
 *    emitter fewer than the best so far. A branch is cut off when the pulses
 *    still missing need more emitters than are left. Branching happens only
 *    on nodes allowed by the minimum cut argument in
 *    PulseNetwork.usefulNextEmitters, so the cut-down choice loses no
 *    solution. When a search ends without finding a smaller set, or the best
 *    set meets the lower bound, the answer is proven optimal.
 *
 * Every candidate set is checked with the maximum flow in PulseNetwork, so any
 * placement this class returns is guaranteed to work.
 *
 * The search is exponential in the worst case, so it runs against a time
 * budget. If the budget runs out the best placement found so far is returned;
 * it is then valid but not proven optimal.
 */
public final class EmitterPlacer {

    /** Candidates evaluated with a full flow update per greedy step. */
    private static final int GREEDY_CANDIDATES = 40;

    private final Instance instance;
    private final Grid grid;
    private final long deadline;
    private final int maxEmitterPulses;
    private boolean timedOut;
    private boolean provenOptimal;
    private List<Integer> forced = new ArrayList<>();
    private int lower;

    /**
     * @param instance     the instance to solve
     * @param budgetMillis time allowed for the branch and bound search
     */
    public EmitterPlacer(Instance instance, long budgetMillis) {
        this.instance = instance;
        this.grid = instance.grid();
        this.deadline = System.currentTimeMillis() + budgetMillis;
        int best = 0;
        for (int v = 0; v < grid.nodeCount(); v++) {
            if (!instance.hasParasite(v)) {
                best = Math.max(best, grid.degree(v));
            }
        }
        this.maxEmitterPulses = best;
    }

    /**
     * Solves the instance.
     *
     * @return the emitter nodes, or null if the instance has no solution
     */
    public List<Integer> solve() {
        if (!isSolvable()) {
            return null;
        }
        forced = forcedEmitters();
        PulseNetwork start = build(forced);
        lower = Math.max(emittersNeeded(instance.totalDemand()),
            forced.size() + emittersNeeded(start.deficit()));

        List<Integer> best = localSearch(greedy(start));

        // Branch and bound, seeded with the greedy result: keep asking for a
        // placement with one emitter fewer than the best so far. When a search
        // finishes without finding one, the best placement is optimal.
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
     * Step 2: the lower bound worked out by solve(). It is the larger of
     * ceil(total demand / most pulses one emitter can send) and the number of
     * forced emitters plus ceil(pulses they cannot deliver / most pulses one
     * emitter can send). The second part is valid because adding an emitter v
     * raises the maximum flow by at most degree(v): removing v's own pulses
     * from any flow leaves a flow for the other emitters.
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
            int before = network.delivered();
            PulseNetwork bestNetwork = null;
            int bestNode = -1;
            for (int v : rankedCandidates(network, null, GREEDY_CANDIDATES)) {
                PulseNetwork trial = network.copy();
                trial.addEmitter(v);
                if (bestNetwork == null || trial.delivered() > bestNetwork.delivered()) {
                    bestNetwork = trial;
                    bestNode = v;
                    if (trial.delivered() - before == grid.degree(v)) {
                        break; // cannot do better than every edge used
                    }
                }
            }
            network = bestNetwork;
            emitters.add(bestNode);
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
        boolean improved = true;
        while (improved) {
            improved = false;
            PulseNetwork network = build(current);
            for (int i = 0; i < current.size(); i++) {
                PulseNetwork trial = network.copy();
                trial.removeEmitter(current.get(i));
                if (trial.isFeasible()) {
                    current.remove(i);
                    improved = true;
                    break;
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
     * Step 4: depth-first branch and bound for a working set of at most k
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
        // Bound: each remaining emitter delivers at most maxEmitterPulses.
        int left = k - chosen.size();
        if (emittersNeeded(network.deficit()) > left) {
            return false;
        }

        // Evaluate every allowed child, then explore the best ones first.
        List<int[]> children = new ArrayList<>();
        for (int v : rankedCandidates(network, forbidden, Integer.MAX_VALUE)) {
            PulseNetwork trial = network.copy();
            trial.addEmitter(v);
            children.add(new int[] {v, trial.delivered()});
        }
        children.sort(Comparator.comparingInt((int[] c) -> -c[1]));

        List<Integer> excluded = new ArrayList<>();
        boolean found = false;
        for (int[] child : children) {
            int v = child[0];
            // A child that cannot reach the target even if every later
            // emitter is perfect is cut off without being expanded.
            int missing = instance.totalDemand() - child[1];
            if (emittersNeeded(missing) > left - 1) {
                continue;
            }
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
     * Free nodes allowed by the minimum cut argument, ordered so the most
     * promising come first: more edges that could carry a new pulse, then
     * closer to a parasite. Returns at most limit nodes.
     */
    private List<Integer> rankedCandidates(PulseNetwork network, boolean[] forbidden, int limit) {
        boolean[] useful = network.usefulNextEmitters();
        boolean[] reach = network.residualReachToSink();
        List<Integer> candidates = new ArrayList<>();
        int n = grid.nodeCount();
        int[] score = new int[n];
        int[] distance = new int[n];
        for (int v = 0; v < n; v++) {
            if (useful[v] && (forbidden == null || !forbidden[v])) {
                candidates.add(v);
                score[v] = network.promisingEdges(v, reach);
                distance[v] = nearestParasite(v);
            }
        }
        Collections.sort(candidates, Comparator
            .comparingInt((Integer v) -> -score[v])
            .thenComparingInt(v -> distance[v])
            .thenComparingInt(v -> v));
        if (candidates.size() > limit) {
            return new ArrayList<>(candidates.subList(0, limit));
        }
        return candidates;
    }

    private int nearestParasite(int v) {
        int best = Integer.MAX_VALUE;
        for (int p : instance.parasites()) {
            best = Math.min(best, grid.distance(v, p));
        }
        return best;
    }
}

package uk.ac.aber.cs31920.assignment;

/**
 * The flow network that checks whether a set of emitters can neutralise every
 * parasite with non-crossing pulses.
 *
 * Pulses must be vertex-disjoint paths (they may only share their emitter or
 * their parasite), which is turned into a flow problem as follows:
 *
 * - Node splitting (22 Oct, slide 190): every grid node v becomes an in-node
 *   and an out-node joined by a "transit" arc of capacity 1, so at most one
 *   pulse can pass through v.
 * - Every grid edge {u, w} becomes the arcs u_out -> w_in and w_out -> u_in,
 *   each of capacity 1, because the grid is undirected (22 Oct, slide 189).
 * - A super source and super sink (22 Oct, slide 191): the source has an arc
 *   to every node's out-node and every node's in-node has an arc to the sink.
 *
 * The capacities depend on what sits on a node:
 *
 * | node     | transit | source -> v_out | v_in -> sink |
 * |----------|---------|-----------------|--------------|
 * | free     | 1       | 0               | 0            |
 * | emitter  | 0       | degree(v)       | 0            |
 * | parasite | 0       | 0               | demand(v)    |
 *
 * So pulses start at emitters, end at parasites and can never pass through
 * either. Since each emitter-to-neighbour arc has capacity 1, an emitter sends
 * at most one pulse per incident edge, exactly as the brief requires. By the
 * integrality of Ford-Fulkerson (22 Oct, slide 195) an integer maximum flow
 * splits into unit paths, which are the pulses. The emitters work if and only
 * if the maximum flow equals the total demand.
 *
 * Emitters can be added and removed while keeping the current flow, so the
 * search only pays for the few augmenting paths that change. A new emitter
 * that no pulse ran through only needs paths that start at it, which are
 * found by a local search around it (FlowNetwork.augmentFrom).
 */
public final class PulseNetwork {

    private final Instance instance;
    private final FlowNetwork network;
    /** transitArc[v], sourceArc[v], sinkArc[v]: the forward arc numbers for node v. */
    private final int[] transitArc;
    private final int[] sourceArc;
    private final int[] sinkArc;
    /** isEmitter[v] is true when v holds an emitter in this network. */
    private final boolean[] isEmitter;
    private int emitterCount;

    /** Builds the network for the instance with no emitters placed yet. */
    public PulseNetwork(Instance instance) {
        this.instance = instance;
        Grid grid = instance.grid();
        int n = grid.nodeCount();
        int source = 2 * n;
        int sink = 2 * n + 1;

        // One transit, source and sink arc per node plus two arcs per edge
        // (the adjacency lists hold every edge twice, once from each side).
        int arcs = 3 * n;
        for (int v = 0; v < n; v++) {
            arcs += grid.degree(v);
        }
        network = new FlowNetwork(2 * n + 2, arcs, source, sink);

        transitArc = new int[n];
        sourceArc = new int[n];
        sinkArc = new int[n];
        isEmitter = new boolean[n];
        for (int v = 0; v < n; v++) {
            boolean parasite = instance.hasParasite(v);
            transitArc[v] = network.addArc(in(v), out(v), parasite ? 0 : 1);
            sourceArc[v] = network.addArc(source, out(v), 0);
            sinkArc[v] = network.addArc(in(v), sink, instance.demand(v));
        }
        for (int v = 0; v < n; v++) {
            for (int w : grid.neighbours(v)) {
                network.addArc(out(v), in(w), 1);
            }
        }
    }

    /** Copy constructor used by the branch and bound search. */
    private PulseNetwork(PulseNetwork other) {
        instance = other.instance;
        network = other.network.copy();
        transitArc = other.transitArc;
        sourceArc = other.sourceArc;
        sinkArc = other.sinkArc;
        isEmitter = other.isEmitter.clone();
        emitterCount = other.emitterCount;
    }

    public PulseNetwork copy() {
        return new PulseNetwork(this);
    }

    private static int in(int v) {
        return 2 * v;
    }

    private static int out(int v) {
        return 2 * v + 1;
    }

    /**
     * Places an emitter on the free node v and recomputes the maximum flow.
     * If a pulse currently runs through v it is cancelled first, because a
     * pulse cannot pass through an emitter; Ford-Fulkerson then reroutes it if
     * it can.
     *
     * @return the number of pulses delivered afterwards
     */
    public int addEmitter(int v) {
        if (isEmitter[v] || instance.hasParasite(v)) {
            throw new IllegalArgumentException("node is not free");
        }
        boolean rerouted = network.flow(transitArc[v]) > 0;
        if (rerouted) {
            network.cancelOneUnit(transitArc[v]);
        }
        network.setCapacity(transitArc[v], 0);
        network.setCapacity(sourceArc[v], instance.grid().degree(v));
        isEmitter[v] = true;
        emitterCount++;
        if (rerouted) {
            // The cancelled pulse freed capacity elsewhere, so new paths may
            // start at other emitters too: run the full algorithm.
            return network.maximise();
        }
        // Otherwise every new pulse must come from v (see augmentFrom).
        return network.augmentFrom(sourceArc[v]);
    }

    /**
     * Exactly how many extra pulses an emitter on the free node v would
     * deliver, without changing this network: the emitter is added as in
     * addEmitter with every change recorded in an undo log, and then the log
     * is rolled back. This costs only the work of the augmenting searches,
     * not a copy of the network.
     */
    public int exactGain(int v) {
        int before = network.value();
        network.beginTrial();
        boolean rerouted = network.flow(transitArc[v]) > 0;
        if (rerouted) {
            network.cancelOneUnit(transitArc[v]);
        }
        network.setCapacity(transitArc[v], 0);
        network.setCapacity(sourceArc[v], instance.grid().degree(v));
        int after = rerouted ? network.maximise() : network.augmentFrom(sourceArc[v]);
        network.rollback();
        return after - before;
    }

    /**
     * Estimates how many extra pulses an emitter on the free node v would
     * deliver, without changing this network. The emitter is tried out with
     * an undo log: a pulse running through v is cancelled, then only v and
     * the emitter whose pulse was cancelled look for new paths, and only
     * through network nodes marked in reach (the nodes that could reach the
     * sink before; see residualReachToSink) or on the cancelled pulse's
     * route. All changes are then rolled back, and reach is left as it was.
     *
     * The estimate is a valid flow, so it never overstates the gain, but it
     * can understate it. It is used only to rank candidates in the greedy
     * construction, never to prove anything. Its cost depends on the size of
     * the reach region, not on the size of the grid.
     */
    public int estimateGain(int v, boolean[] reach) {
        int before = network.value();
        network.beginTrial();
        int freed = -1;
        int[] opened = new int[0];
        if (network.flow(transitArc[v]) > 0) {
            freed = network.cancelOneUnit(transitArc[v]);
            // The cancelled pulse's route is free again, so allow it too.
            opened = network.markTrialNodes(reach);
        }
        network.setCapacity(transitArc[v], 0);
        network.setCapacity(sourceArc[v], instance.grid().degree(v));
        network.augmentFrom(sourceArc[v], reach);
        if (freed >= 0) {
            network.augmentFrom(freed, reach);
        }
        int gain = network.value() - before;
        network.rollback();
        for (int w : opened) {
            reach[w] = false;
        }
        return gain;
    }

    /**
     * Removes the emitter on v, cancelling its pulses, and recomputes the
     * maximum flow. Only meant for a network whose flow was maximum before.
     *
     * @return the number of pulses delivered afterwards
     */
    public int removeEmitter(int v) {
        if (!isEmitter[v]) {
            throw new IllegalArgumentException("no emitter on node");
        }
        while (network.flow(sourceArc[v]) > 0) {
            network.cancelOneUnit(sourceArc[v]);
        }
        network.setCapacity(sourceArc[v], 0);
        network.setCapacity(transitArc[v], 1);
        isEmitter[v] = false;
        emitterCount--;
        // Every augmenting path now has to end at a parasite that lost a
        // pulse, since all others are still full. So search backwards from
        // those parasites until a whole round finds nothing.
        boolean progress = true;
        while (progress) {
            progress = false;
            for (int p : instance.parasites()) {
                if (network.flow(sinkArc[p]) < instance.demand(p) && network.augmentInto(sinkArc[p])) {
                    progress = true;
                }
            }
        }
        return network.value();
    }

    // ----- Relaxed networks, used only for lower bounds (see RegionBound) -----

    /**
     * Lets the free node v act as a source of up to degree(v) pulses while it
     * also stays available as a transit node. A network in which several
     * nodes are opened like this is a relaxation: every real placement that
     * only uses opened nodes as emitters has a flow that is also valid here.
     * Call maximiseRelaxed afterwards (or augmentRelaxed for one node when
     * the flow was already maximum).
     */
    public void openRelaxedSource(int v) {
        network.setCapacity(sourceArc[v], instance.grid().degree(v));
    }

    /** Maximum flow of a relaxed network from scratch (Dinitz). */
    public int maximiseRelaxed() {
        return network.maximise();
    }

    /** Opens v and updates a maximum flow with paths that start at v only. */
    public int augmentRelaxed(int v) {
        openRelaxedSource(v);
        return network.augmentFrom(sourceArc[v]);
    }

    /** Starts recording changes so that rollbackTrial can undo them. */
    public void beginTrial() {
        network.beginTrial();
    }

    /** Keeps every change since beginTrial. */
    public void commitTrial() {
        network.commit();
    }

    /** Undoes every change since beginTrial. */
    public void rollbackTrial() {
        network.rollback();
    }

    /**
     * Pulses the emitters other than v could still send: the most that
     * removing v could ever be made up for.
     */
    public int spareCapacityWithout(int v) {
        int spare = 0;
        Grid grid = instance.grid();
        for (int w = 0; w < grid.nodeCount(); w++) {
            if (isEmitter[w] && w != v) {
                spare += grid.degree(w) - network.flow(sourceArc[w]);
            }
        }
        return spare;
    }

    /** Pulses sent by the emitter on v in the current flow. */
    public int pulsesFrom(int v) {
        return network.flow(sourceArc[v]);
    }

    /** The number of pulses the current emitters deliver (the maximum flow). */
    public int delivered() {
        return network.value();
    }

    /** Pulses still missing to neutralise every parasite. */
    public int deficit() {
        return instance.totalDemand() - network.value();
    }

    public boolean isFeasible() {
        return deficit() == 0;
    }

    public boolean isEmitter(int v) {
        return isEmitter[v];
    }

    public int emitterCount() {
        return emitterCount;
    }

    /** True if v holds neither a parasite nor an emitter. */
    public boolean isFree(int v) {
        return !isEmitter[v] && !instance.hasParasite(v);
    }

    /**
     * Marks the free nodes that are worth trying as the next emitter.
     *
     * After the maximum flow, let Z be the set of network nodes that can still
     * reach the sink in the residual network. The arcs into Z form a minimum
     * cut whose capacity is the current flow (max-flow min-cut, 23 Oct,
     * slides 215-219). A new emitter at w only adds an arc source -> w_out, and
     * that arc crosses the cut only if w_out is in Z. Removing transit arcs can
     * only lower the cut. So if every new emitter had w_out outside Z the cut,
     * and with it the flow, could not grow. Every set of extra emitters that
     * fixes the deficit must therefore contain at least one node marked here.
     * The branch and bound search uses this to branch only on these nodes.
     *
     * @param reach the result of residualReachToSink() for the current flow
     */
    public boolean[] usefulNextEmitters(boolean[] reach) {
        int n = instance.grid().nodeCount();
        boolean[] useful = new boolean[n];
        for (int v = 0; v < n; v++) {
            useful[v] = isFree(v) && reach[out(v)];
        }
        return useful;
    }

    /**
     * Number of neighbours of v from which an extra pulse could still go on to
     * a parasite. Used only to order candidates (a heuristic, not a bound).
     */
    public int promisingEdges(int v, boolean[] reach) {
        int count = 0;
        for (int w : instance.grid().neighbours(v)) {
            if (!isEmitter[w] && reach[in(w)]) {
                count++;
            }
        }
        return count;
    }

    /** Residual reachability to the sink, indexed by network node. */
    public boolean[] residualReachToSink() {
        return network.canReachSink();
    }
}

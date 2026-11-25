package uk.ac.aber.cs31920.assignment;

import java.util.ArrayList;
import java.util.List;

/**
 * The flow network that checks whether a set of emitters can neutralise every
 * parasite with non-crossing pulses.
 *
 * Pulses must be paths that share no node (they may only share their emitter
 * or their parasite), which is turned into a flow problem as follows:
 *
 * - Node capacities (22 Oct, slide 190): every grid node v is replaced by an
 *   in-node and an out-node joined by a "transit" edge of capacity 1, so at most
 *   one pulse can pass through v.
 * - Every undirected grid edge {u, w} is replaced by the directed edges
 *   u_out -> w_in and w_out -> u_in, each of capacity 1 (undirected edges,
 *   22 Oct, slide 189).
 * - A super source and super sink (22 Oct, slide 191): the super source has an
 *   edge to every node's out-node and every node's in-node has an edge to the
 *   super sink.
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
 * either. Since each emitter-to-neighbour edge has capacity 1, an emitter sends
 * at most one pulse per incident edge, exactly as the brief requires. Since all
 * capacities are integers, the Ford/Fulkerson algorithm gives an integer
 * maximum flow (22 Oct, slide 195), which splits into s-t paths of one
 * unit each: the pulses. The emitters work if and only if the value of the
 * maximum flow equals the total demand.
 *
 * Emitters can be added and removed while keeping the current flow, so the
 * search only pays for the few augmenting paths that change. A new emitter
 * that no pulse ran through only needs paths that start at it, which are
 * found by a local search around it (FlowNetwork.augmentFrom).
 */
public final class PulseNetwork {

    private final Instance instance;
    private final FlowNetwork network;
    /** transitEdge[v], sourceEdge[v], sinkEdge[v]: the forward edge numbers for node v. */
    private final int[] transitEdge;
    private final int[] sourceEdge;
    private final int[] sinkEdge;
    /** isEmitter[v] is true when v holds an emitter in this network. */
    private final boolean[] isEmitter;

    /** Builds the network for the instance with no emitters placed yet. */
    public PulseNetwork(Instance instance) {
        this.instance = instance;
        Grid grid = instance.grid();
        int n = grid.nodeCount();
        int source = 2 * n;
        int sink = 2 * n + 1;

        // One transit, source and sink edge per node plus two directed edges per
        // grid edge (the adjacency lists hold every grid edge twice, once from each side).
        int edges = 3 * n;
        for (int v = 0; v < n; v++) {
            edges += grid.degree(v);
        }
        network = new FlowNetwork(2 * n + 2, edges, source, sink);

        transitEdge = new int[n];
        sourceEdge = new int[n];
        sinkEdge = new int[n];
        isEmitter = new boolean[n];
        for (int v = 0; v < n; v++) {
            boolean parasite = instance.hasParasite(v);
            transitEdge[v] = network.addEdge(in(v), out(v), parasite ? 0 : 1);
            sourceEdge[v] = network.addEdge(source, out(v), 0);
            sinkEdge[v] = network.addEdge(in(v), sink, instance.demand(v));
        }
        for (int v = 0; v < n; v++) {
            for (int w : grid.neighbours(v)) {
                network.addEdge(out(v), in(w), 1);
            }
        }
    }

    /** Copy constructor used by the branch and bound search. */
    private PulseNetwork(PulseNetwork other) {
        instance = other.instance;
        network = other.network.copy();
        transitEdge = other.transitEdge;
        sourceEdge = other.sourceEdge;
        sinkEdge = other.sinkEdge;
        isEmitter = other.isEmitter.clone();
    }

    /** An independent copy with the same emitters and flow. */
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
     * pulse cannot pass through an emitter; the Ford/Fulkerson algorithm then reroutes it if
     * it can.
     *
     * @return the number of pulses delivered afterwards
     */
    public int addEmitter(int v) {
        if (!isFree(v)) {
            throw new IllegalArgumentException("node is not free");
        }
        isEmitter[v] = true;
        return placeAndMaximise(v);
    }

    /** Turns v into a source in the flow network and restores a maximum flow. */
    private int placeAndMaximise(int v) {
        boolean rerouted = network.flow(transitEdge[v]) > 0;
        if (rerouted) {
            network.cancelOneUnit(transitEdge[v]);
        }
        network.setCapacity(transitEdge[v], 0);
        network.setCapacity(sourceEdge[v], instance.grid().degree(v));
        if (rerouted) {
            // The cancelled pulse freed capacity elsewhere, so new paths may
            // start at other emitters too: run the full algorithm.
            return network.maximise();
        }
        // Otherwise every new pulse must come from v (see augmentFrom).
        return network.augmentFrom(sourceEdge[v]);
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
        int after = placeAndMaximise(v);
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
     * can understate it. It is used only to rank candidates in the constructive
     * heuristic, never to prove anything. Its cost depends on the size of
     * the reach region, not on the size of the grid.
     */
    public int estimateGain(int v, boolean[] reach) {
        int before = network.value();
        network.beginTrial();
        int freed = -1;
        int[] opened = new int[0];
        if (network.flow(transitEdge[v]) > 0) {
            freed = network.cancelOneUnit(transitEdge[v]);
            // The cancelled pulse's route is free again, so allow it too.
            opened = network.markTrialNodes(reach);
        }
        network.setCapacity(transitEdge[v], 0);
        network.setCapacity(sourceEdge[v], instance.grid().degree(v));
        network.augmentFrom(sourceEdge[v], reach);
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
        while (network.flow(sourceEdge[v]) > 0) {
            network.cancelOneUnit(sourceEdge[v]);
        }
        network.setCapacity(sourceEdge[v], 0);
        network.setCapacity(transitEdge[v], 1);
        isEmitter[v] = false;
        // Every augmenting path now has to end at a parasite that lost a
        // pulse, since all others are still full. So search backwards from
        // those parasites until a whole round finds nothing.
        boolean progress = true;
        while (progress) {
            progress = false;
            for (int p : instance.parasites()) {
                if (network.flow(sinkEdge[p]) < instance.demand(p) && network.augmentInto(sinkEdge[p])) {
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
        network.setCapacity(sourceEdge[v], instance.grid().degree(v));
    }

    /** Maximum flow of a relaxed network from scratch (Dinitz). */
    public int maximiseRelaxed() {
        return network.maximise();
    }

    /** Opens v and updates a maximum flow with paths that start at v only. */
    public int augmentRelaxed(int v) {
        openRelaxedSource(v);
        return network.augmentFrom(sourceEdge[v]);
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
                spare += grid.degree(w) - network.flow(sourceEdge[w]);
            }
        }
        return spare;
    }

    /** Pulses sent by the emitter on v in the current flow. */
    public int pulsesFrom(int v) {
        return network.flow(sourceEdge[v]);
    }

    /** Pulses still missing to neutralise every parasite. */
    public int deficit() {
        return instance.totalDemand() - network.value();
    }

    /** True if the emitters neutralise every parasite. */
    public boolean isFeasible() {
        return deficit() == 0;
    }

    /** True if v holds neither a parasite nor an emitter. */
    public boolean isFree(int v) {
        return !isEmitter[v] && !instance.hasParasite(v);
    }

    /**
     * Marks the free nodes that are worth trying as the next emitter.
     *
     * After the maximum flow, let Z be the set of network nodes that can still
     * reach the sink in the residual network. The edges into Z form a minimum
     * cut whose capacity is the current flow (Max-Flow Min-Cut Theorem, 23 Oct,
     * slides 215-219). A new emitter at w only adds an edge source -> w_out, and
     * that edge crosses the cut only if w_out is in Z. Removing transit edges can
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

    /** Pulses the parasite on p still lacks in the current flow. */
    public int shortfall(int p) {
        return instance.demand(p) - network.flow(sinkEdge[p]);
    }

    /**
     * The free nodes that could help the parasite on p, which must be short
     * of pulses: the nodes w whose out-node can reach p's in-node in the
     * residual network without passing through the sink.
     *
     * Every set of extra emitters that neutralises all parasites contains
     * one of these nodes. Let Y be the network nodes that can reach p_in in
     * this way. No residual edge enters Y from outside except from the source
     * or sink, so every edge into Y (other than from the source) is full and
     * every edge out of Y (other than into the sink) is empty. The pulses that
     * the parasites in Y receive therefore equal the capacity into Y. Adding
     * emitters outside Y can only remove capacity into Y (a transit edge), so
     * those parasites could never receive more pulses in total, yet p needs
     * more. This is the Max-Flow Min-Cut Theorem (23 Oct, slides 215-219)
     * applied to the cut around Y.
     */
    public List<Integer> helpersOf(int p) {
        boolean[] reach = network.canReach(in(p));
        List<Integer> helpers = new ArrayList<>();
        for (int w = 0; w < instance.grid().nodeCount(); w++) {
            if (isFree(w) && reach[out(w)]) {
                helpers.add(w);
            }
        }
        return helpers;
    }

    /** Residual reachability to the sink, indexed by network node. */
    public boolean[] residualReachToSink() {
        return network.canReachSink();
    }
}

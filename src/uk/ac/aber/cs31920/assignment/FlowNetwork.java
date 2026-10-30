package uk.ac.aber.cs31920.assignment;

import java.util.Arrays;

/**
 * A directed flow network with integer capacities and the Ford-Fulkerson
 * method for maximum flow (lecture 22 Oct, slide 203; runtime 23 Oct, slide 212).
 *
 * Augmenting paths are found by breadth-first search in the residual network,
 * so each augmentation costs O(|V| + |E|). Flow values here are at most the
 * total parasite demand, which keeps the number of augmentations small.
 *
 * Every arc a is stored together with its reverse arc a ^ 1 (forward arcs have
 * even numbers). The reverse arc has capacity 0 and always carries minus the
 * flow of its partner, so the residual capacity of any arc is cap - flow.
 *
 * The topology (head, next, to) never changes after construction, so copies
 * share it and only clone the capacity and flow arrays. That makes copying a
 * network cheap, which the branch and bound search relies on.
 */
public final class FlowNetwork {

    private final int nodeCount;
    private final int source;
    private final int sink;

    // Topology, shared between copies.
    private int[] head;
    private int[] next;
    private int[] to;
    private int arcCount;

    // State, cloned by copies.
    private int[] cap;
    private int[] flow;
    private int value;

    // Scratch space for the breadth-first search, shared between copies (the
    // search is single threaded, so only one network searches at a time).
    private int[] parentArc;
    private int[] queue;
    private int[] visitMark;
    private int[] visitStamp;

    /**
     * Creates an empty network.
     *
     * @param nodeCount number of nodes, numbered 0 .. nodeCount-1
     * @param maxArcs   number of forward arcs that will be added
     * @param source    the source node
     * @param sink      the sink node
     */
    public FlowNetwork(int nodeCount, int maxArcs, int source, int sink) {
        this.nodeCount = nodeCount;
        this.source = source;
        this.sink = sink;
        head = new int[nodeCount];
        Arrays.fill(head, -1);
        next = new int[2 * maxArcs];
        to = new int[2 * maxArcs];
        cap = new int[2 * maxArcs];
        flow = new int[2 * maxArcs];
        parentArc = new int[nodeCount];
        queue = new int[nodeCount];
        visitMark = new int[nodeCount];
        visitStamp = new int[1];
    }

    /** Copy constructor: shares topology and scratch space, clones capacities and flow. */
    private FlowNetwork(FlowNetwork other) {
        nodeCount = other.nodeCount;
        source = other.source;
        sink = other.sink;
        head = other.head;
        next = other.next;
        to = other.to;
        arcCount = other.arcCount;
        cap = other.cap.clone();
        flow = other.flow.clone();
        value = other.value;
        parentArc = other.parentArc;
        queue = other.queue;
        visitMark = other.visitMark;
        visitStamp = other.visitStamp;
    }

    public FlowNetwork copy() {
        return new FlowNetwork(this);
    }

    /**
     * Adds an arc u -> v with the given capacity (and its reverse arc).
     *
     * @return the number of the forward arc
     */
    public int addArc(int u, int v, int capacity) {
        int a = 2 * arcCount++;
        link(a, u, v, capacity);
        link(a + 1, v, u, 0);
        return a;
    }

    private void link(int a, int u, int v, int capacity) {
        to[a] = v;
        cap[a] = capacity;
        next[a] = head[u];
        head[u] = a;
    }

    public int value() {
        return value;
    }

    public int flow(int arc) {
        return flow[arc];
    }

    /**
     * Changes the capacity of a forward arc. The caller must make sure the arc
     * does not carry more flow than the new capacity (see cancelOneUnit).
     */
    public void setCapacity(int arc, int capacity) {
        if (flow[arc] > capacity) {
            throw new IllegalStateException("arc carries more flow than the new capacity");
        }
        cap[arc] = capacity;
    }

    /**
     * Runs Ford-Fulkerson from the current flow until no augmenting path is
     * left. Because the method can start from any valid flow, a network that
     * gained an emitter only has to find the few extra paths.
     *
     * @return the maximum flow value
     */
    public int maximise() {
        while (augment()) {
            // keep augmenting
        }
        return value;
    }

    /**
     * Finds one augmenting path with breadth-first search and pushes the
     * bottleneck amount of flow along it.
     *
     * @return false if the sink is not reachable in the residual network
     */
    private boolean augment() {
        int stamp = ++visitStamp[0];
        int front = 0;
        int back = 0;
        queue[back++] = source;
        visitMark[source] = stamp;
        boolean found = false;
        while (front < back && !found) {
            int u = queue[front++];
            for (int a = head[u]; a >= 0; a = next[a]) {
                int v = to[a];
                if (visitMark[v] != stamp && cap[a] - flow[a] > 0) {
                    visitMark[v] = stamp;
                    parentArc[v] = a;
                    if (v == sink) {
                        found = true;
                        break;
                    }
                    queue[back++] = v;
                }
            }
        }
        if (!found) {
            return false;
        }

        // Walk back from the sink to find the bottleneck, then push.
        int bottleneck = Integer.MAX_VALUE;
        for (int v = sink; v != source; v = to[parentArc[v] ^ 1]) {
            int a = parentArc[v];
            bottleneck = Math.min(bottleneck, cap[a] - flow[a]);
        }
        for (int v = sink; v != source; v = to[parentArc[v] ^ 1]) {
            int a = parentArc[v];
            flow[a] += bottleneck;
            flow[a ^ 1] -= bottleneck;
        }
        value += bottleneck;
        return true;
    }

    /**
     * Removes one unit of flow that passes through the given forward arc,
     * keeping flow conservation everywhere else. The unit is traced forwards
     * to the sink and backwards to the source (or round a cycle back to the
     * arc) and every arc on the way loses one unit.
     *
     * This is how the network reacts when a node that carries a pulse becomes
     * an emitter, or an emitter is removed: the affected pulse is cancelled
     * and Ford-Fulkerson is free to find a new route afterwards.
     */
    public void cancelOneUnit(int arc) {
        if (flow[arc] <= 0) {
            throw new IllegalStateException("arc carries no flow");
        }
        push(arc, -1);
        int tail = to[arc ^ 1];

        // Forwards: the head now receives one unit less than it sends.
        int v = to[arc];
        boolean closedCycle = false;
        while (v != sink) {
            if (v == tail) {
                closedCycle = true;
                break;
            }
            int out = outgoingFlowArc(v);
            push(out, -1);
            v = to[out];
        }

        // Backwards: the tail now sends one unit less than it receives.
        if (!closedCycle) {
            int u = tail;
            while (u != source) {
                int in = incomingFlowArc(u);
                push(in, -1);
                u = to[in ^ 1];
            }
            value -= 1;
        }
    }

    private void push(int arc, int amount) {
        flow[arc] += amount;
        flow[arc ^ 1] -= amount;
    }

    /** A forward arc leaving u that carries positive flow. */
    private int outgoingFlowArc(int u) {
        for (int a = head[u]; a >= 0; a = next[a]) {
            if ((a & 1) == 0 && flow[a] > 0) {
                return a;
            }
        }
        throw new IllegalStateException("flow conservation violated");
    }

    /** A forward arc entering u that carries positive flow. */
    private int incomingFlowArc(int u) {
        for (int a = head[u]; a >= 0; a = next[a]) {
            if ((a & 1) == 1 && flow[a ^ 1] > 0) {
                return a ^ 1;
            }
        }
        throw new IllegalStateException("flow conservation violated");
    }

    /**
     * Marks every node that can still reach the sink in the residual network.
     * After a maximum flow this is the sink side of a minimum cut
     * (max-flow min-cut theorem, 23 Oct, slides 215-219).
     */
    public boolean[] canReachSink() {
        boolean[] reach = new boolean[nodeCount];
        int front = 0;
        int back = 0;
        queue[back++] = sink;
        reach[sink] = true;
        while (front < back) {
            int v = queue[front++];
            // An arc u -> v is usable if it has residual capacity; its reverse
            // arc (stored at v) points back to u.
            for (int r = head[v]; r >= 0; r = next[r]) {
                int a = r ^ 1;
                int u = to[r];
                if (!reach[u] && cap[a] - flow[a] > 0) {
                    reach[u] = true;
                    queue[back++] = u;
                }
            }
        }
        return reach;
    }
}

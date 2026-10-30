package uk.ac.aber.cs31920.assignment;

import java.util.Arrays;

/**
 * A directed flow network with integer capacities and Dinitz's algorithm for
 * maximum flow (lecture 28 Oct, slides 229-234).
 *
 * Each phase builds the layered residual network with one breadth-first
 * search and then adds a blocking flow, found by depth-first search with a
 * "current arc" pointer per node so that no arc is looked at twice after it
 * has turned out to be useless in that phase (slide 232).
 *
 * Why Dinitz rather than plain Ford-Fulkerson (22 Oct, slide 203)? Ford-
 * Fulkerson needs up to |f| rounds of O(n + e) each (23 Oct, slide 212), and
 * here |f| is the total demand, which can grow linearly with the grid. The
 * general bound for Dinitz, O(n^2 e) (slide 234), looks worse, but the pulse
 * network is almost a unit network: after node splitting a free node carries
 * at most one unit, and an emitter or parasite at most 8 (its degree). So the
 * remaining flow splits into paths that share each node at most 8 times.
 * After d phases every augmenting path is longer than d, so at most 8n/d
 * units of flow are still missing. Taking d = sqrt(n) shows there are only
 * O(sqrt(n)) phases. A blocking flow takes O(e) time in such a network,
 * because a path saturates its unit arcs and each arc is then dropped. The
 * total
 * is O(sqrt(n) * e) instead of O(|f| * e).
 *
 * Every arc a is stored together with its reverse arc a ^ 1 (forward arcs have
 * even numbers). The reverse arc has capacity 0 and always carries minus the
 * flow of its partner, so the residual capacity of any arc is cap - flow.
 *
 * The topology (head, next, to) never changes after construction, so copies
 * share it and only clone the capacity and flow arrays. That makes copying a
 * network cheap, which the branch and bound search relies on. Any valid flow
 * can be the starting point, so after a small change (one more emitter) only
 * the few missing paths have to be found.
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

    // Scratch space for the searches, shared between copies (the program is
    // single threaded, so only one network searches at a time).
    private int[] level;
    private int[] currentArc;
    private int[] pathArc;
    private int[] queue;
    private int[] parentArc;
    private int[] visitMark;
    private int[] visitStamp;

    // Undo log for trial changes (see beginTrial / rollback). An entry
    // (arc, amount) undoes a push; an entry (~arc, old) undoes a capacity change.
    private int[] trailArc = new int[64];
    private int[] trailValue = new int[64];
    private int trailSize = -1;
    private int valueBeforeTrial;

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
        level = new int[nodeCount];
        currentArc = new int[nodeCount];
        pathArc = new int[nodeCount];
        queue = new int[nodeCount];
        parentArc = new int[nodeCount];
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
        level = other.level;
        currentArc = other.currentArc;
        pathArc = other.pathArc;
        queue = other.queue;
        parentArc = other.parentArc;
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
        record(~arc, cap[arc]);
        cap[arc] = capacity;
    }

    /**
     * Runs Dinitz's algorithm from the current flow: while the sink can be
     * reached in the residual network, build the layered network and add a
     * blocking flow (28 Oct, slide 231).
     *
     * @return the maximum flow value
     */
    public int maximise() {
        while (buildLevels()) {
            System.arraycopy(head, 0, currentArc, 0, nodeCount);
            int pushed;
            do {
                pushed = pushPathInLevels();
                value += pushed;
            } while (pushed > 0);
        }
        return value;
    }

    /**
     * Breadth-first search from the source that gives every node its distance
     * in the residual network. Nodes at the sink's distance or beyond are not
     * expanded, since no shortest path uses them.
     *
     * @return false if the sink cannot be reached
     */
    private boolean buildLevels() {
        Arrays.fill(level, -1);
        int front = 0;
        int back = 0;
        queue[back++] = source;
        level[source] = 0;
        while (front < back) {
            int u = queue[front++];
            if (level[sink] >= 0 && level[u] >= level[sink]) {
                break;
            }
            for (int a = head[u]; a >= 0; a = next[a]) {
                int v = to[a];
                if (level[v] < 0 && cap[a] - flow[a] > 0) {
                    level[v] = level[u] + 1;
                    queue[back++] = v;
                }
            }
        }
        return level[sink] >= 0;
    }

    /**
     * Finds one source-sink path in the layered network and pushes its
     * bottleneck along it. The search is iterative (paths in a large grid can
     * be long) and advances each node's current arc past arcs that are full
     * or lead to a dead end, so over one phase every arc is skipped at most
     * once (slide 232). A dead-end node is removed from the layers.
     *
     * @return the amount pushed, or 0 if the flow in the layers is blocking
     */
    private int pushPathInLevels() {
        int depth = 0;
        int u = source;
        while (true) {
            if (u == sink) {
                int bottleneck = Integer.MAX_VALUE;
                for (int i = 0; i < depth; i++) {
                    int a = pathArc[i];
                    bottleneck = Math.min(bottleneck, cap[a] - flow[a]);
                }
                for (int i = 0; i < depth; i++) {
                    push(pathArc[i], bottleneck);
                }
                return bottleneck;
            }
            // Advance u's current arc to the next usable arc of the layers.
            int a = currentArc[u];
            while (a >= 0 && (level[to[a]] != level[u] + 1 || cap[a] - flow[a] <= 0)) {
                a = next[a];
            }
            currentArc[u] = a;
            if (a >= 0) {
                pathArc[depth++] = a;
                u = to[a];
            } else {
                // Dead end: drop u from the layers and step back.
                if (u == source) {
                    return 0;
                }
                level[u] = -1;
                int back = pathArc[--depth];
                u = to[back ^ 1];
                currentArc[u] = next[currentArc[u]];
            }
        }
    }

    /**
     * Ford-Fulkerson (22 Oct, slide 203) restricted to augmenting paths that
     * start with the given arc out of the source.
     *
     * This is enough after one source arc has been added to a network whose
     * flow was already maximum. Let X be the nodes the source could reach
     * before. No residual arc left X, otherwise the flow would not have been
     * maximum. A simple augmenting path that entered X could therefore never
     * leave it again, so every augmenting path, now and after each later
     * augmentation, starts with the new arc and stays outside X. Searching
     * from the head of that arc finds them all, and the breadth-first search
     * stops as soon as it meets the sink, so it usually only looks at the
     * part of the grid around the new emitter instead of the whole network.
     *
     * @param startArc a forward arc leaving the source
     * @return the maximum flow value
     */
    public int augmentFrom(int startArc) {
        return augmentFrom(startArc, null);
    }

    /**
     * As augmentFrom(startArc), but the search may only visit nodes marked in
     * allowed (null allows all). With a restriction the result is a valid flow
     * that may be below the maximum; it is used for quick estimates.
     */
    public int augmentFrom(int startArc, boolean[] allowed) {
        int start = to[startArc];
        while (cap[startArc] - flow[startArc] > 0) {
            int stamp = ++visitStamp[0];
            int front = 0;
            int back = 0;
            queue[back++] = start;
            visitMark[start] = stamp;
            visitMark[source] = stamp; // never walk back into the source
            boolean found = false;
            while (front < back && !found) {
                int u = queue[front++];
                for (int a = head[u]; a >= 0; a = next[a]) {
                    int v = to[a];
                    if (visitMark[v] != stamp && cap[a] - flow[a] > 0
                            && (allowed == null || allowed[v])) {
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
                break;
            }
            // Every arc between grid nodes has capacity 1, so push one unit.
            for (int v = sink; v != start; v = to[parentArc[v] ^ 1]) {
                push(parentArc[v], 1);
            }
            push(startArc, 1);
            value += 1;
        }
        return value;
    }

    /**
     * The mirror image of augmentFrom: Ford-Fulkerson restricted to
     * augmenting paths that end with the given arc into the sink. The search
     * runs backwards from the tail of that arc to the source and stops when
     * it gets there.
     *
     * @param endArc a forward arc entering the sink
     * @return true if at least one unit was pushed
     */
    public boolean augmentInto(int endArc) {
        int end = to[endArc ^ 1];
        boolean pushedAny = false;
        while (cap[endArc] - flow[endArc] > 0) {
            int stamp = ++visitStamp[0];
            int front = 0;
            int back = 0;
            queue[back++] = end;
            visitMark[end] = stamp;
            visitMark[sink] = stamp; // never walk back into the sink
            boolean found = false;
            while (front < back && !found) {
                int v = queue[front++];
                // Arc u -> v is usable if it has residual capacity; it is the
                // partner of the arc r stored at v.
                for (int r = head[v]; r >= 0; r = next[r]) {
                    int a = r ^ 1;
                    int u = to[r];
                    if (visitMark[u] != stamp && cap[a] - flow[a] > 0) {
                        visitMark[u] = stamp;
                        parentArc[u] = a; // the arc leaving u towards the end
                        if (u == source) {
                            found = true;
                            break;
                        }
                        queue[back++] = u;
                    }
                }
            }
            if (!found) {
                break;
            }
            int bottleneck = cap[endArc] - flow[endArc];
            for (int u = source; u != end; u = to[parentArc[u]]) {
                int a = parentArc[u];
                bottleneck = Math.min(bottleneck, cap[a] - flow[a]);
            }
            for (int u = source; u != end; u = to[parentArc[u]]) {
                push(parentArc[u], bottleneck);
            }
            push(endArc, bottleneck);
            value += bottleneck;
            pushedAny = true;
        }
        return pushedAny;
    }

    /**
     * Removes one unit of flow that passes through the given forward arc,
     * keeping flow conservation everywhere else. The unit is traced forwards
     * to the sink and backwards to the source (or round a cycle back to the
     * arc) and every arc on the way loses one unit.
     *
     * This is how the network reacts when a node that carries a pulse becomes
     * an emitter, or an emitter is removed: the affected pulse is cancelled
     * and the next maximise call is free to find a new route.
     *
     * @return the source arc the cancelled unit started from, or -1 if the
     *         unit went round a cycle
     */
    public int cancelOneUnit(int arc) {
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
        if (closedCycle) {
            return -1;
        }
        int startArc = arc;
        int u = tail;
        while (u != source) {
            startArc = incomingFlowArc(u);
            push(startArc, -1);
            u = to[startArc ^ 1];
        }
        value -= 1;
        return startArc;
    }

    private void push(int arc, int amount) {
        record(arc, amount);
        flow[arc] += amount;
        flow[arc ^ 1] -= amount;
    }

    /**
     * Starts recording changes so that rollback can undo them. Trying out a
     * candidate emitter then costs only as much as the paths it changes,
     * instead of O(e) for copying the whole network.
     */
    public void beginTrial() {
        trailSize = 0;
        valueBeforeTrial = value;
    }

    /** Keeps every change since beginTrial and stops recording. */
    public void commit() {
        trailSize = -1;
    }

    /** Undoes every change since beginTrial and stops recording. */
    public void rollback() {
        for (int i = trailSize - 1; i >= 0; i--) {
            int arc = trailArc[i];
            if (arc >= 0) {
                flow[arc] -= trailValue[i];
                flow[arc ^ 1] += trailValue[i];
            } else {
                cap[~arc] = trailValue[i];
            }
        }
        trailSize = -1;
        value = valueBeforeTrial;
    }

    /**
     * Adds to mark every node at either end of an arc whose flow changed
     * since beginTrial, and returns the nodes that were newly marked so the
     * caller can clear them again.
     */
    public int[] markTrialNodes(boolean[] mark) {
        int[] added = new int[2 * Math.max(trailSize, 0)];
        int count = 0;
        for (int i = 0; i < trailSize; i++) {
            int arc = trailArc[i];
            if (arc < 0) {
                continue; // a capacity change, not a flow change
            }
            for (int v : new int[] {to[arc], to[arc ^ 1]}) {
                if (!mark[v]) {
                    mark[v] = true;
                    added[count++] = v;
                }
            }
        }
        return Arrays.copyOf(added, count);
    }

    private void record(int arc, int amount) {
        if (trailSize < 0) {
            return; // not in a trial
        }
        if (trailSize == trailArc.length) {
            trailArc = Arrays.copyOf(trailArc, 2 * trailSize);
            trailValue = Arrays.copyOf(trailValue, 2 * trailSize);
        }
        trailArc[trailSize] = arc;
        trailValue[trailSize] = amount;
        trailSize++;
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

package uk.ac.aber.cs31920.assignment;

import java.util.Arrays;

/**
 * A directed flow network with integer capacities and Dinitz's algorithm for
 * maximum flow (lecture 28 Oct, slides 229-234).
 * Each round computes the layered residual network with one breadth first
 * search (BFS) and then adds a blocking flow, found by depth first search
 * (DFS) that removes dead ends as they are discovered (slide 232). A
 * "current edge" pointer per node makes sure no edge is looked at twice
 * after it has turned out to be useless in that round.
 * Why Dinitz's algorithm rather than the Ford/Fulkerson algorithm (22 Oct,
 * slide 203)? The Ford/Fulkerson algorithm has pseudo-polynomial runtime
 * (23 Oct, slide 212): it may need as many rounds as the value of the
 * maximum flow, here the total demand. Dinitz's algorithm has polynomial
 * runtime O(n^2 e) (28 Oct, slide 234), independent of the capacities.
 * Every edge a is stored together with its backward edge a ^ 1 (forward
 * edges have even numbers). The backward edge has capacity 0 and always
 * carries minus the flow of its partner, so the residual capacity of any
 * edge is cap - flow.
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
    private final int[] head;
    private final int[] next;
    private final int[] to;
    private int edgeCount;

    // State, cloned by copies.
    private final int[] cap;
    private final int[] flow;
    private int value;

    // Scratch space for the searches, shared between copies (the program is
    // single threaded, so only one network searches at a time).
    private final int[] distance;
    private final int[] currentEdge;
    private final int[] pathEdge;
    private final int[] queue;
    private final int[] parentEdge;
    private final int[] visitMark;
    /** One counter in an array so that copies share it with visitMark. */
    private final int[] visitStamp;

    // Undo log for trial changes (see beginTrial / rollback). An entry
    // (edge, amount) undoes a push; an entry (~edge, old) undoes a capacity change.
    private int[] trailEdge = new int[64];
    private int[] trailValue = new int[64];
    private int trailSize = -1;
    private int valueBeforeTrial;

    /**
     * Creates an empty network.
     *
     * @param nodeCount number of nodes, numbered 0 .. nodeCount-1
     * @param maxEdges   number of forward edges that will be added
     * @param source    the source node
     * @param sink      the sink node
     */
    public FlowNetwork(int nodeCount, int maxEdges, int source, int sink) {
        this.nodeCount = nodeCount;
        this.source = source;
        this.sink = sink;
        head = new int[nodeCount];
        Arrays.fill(head, -1);
        next = new int[2 * maxEdges];
        to = new int[2 * maxEdges];
        cap = new int[2 * maxEdges];
        flow = new int[2 * maxEdges];
        distance = new int[nodeCount];
        currentEdge = new int[nodeCount];
        pathEdge = new int[nodeCount];
        queue = new int[nodeCount];
        parentEdge = new int[nodeCount];
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
        edgeCount = other.edgeCount;
        cap = other.cap.clone();
        flow = other.flow.clone();
        value = other.value;
        distance = other.distance;
        currentEdge = other.currentEdge;
        pathEdge = other.pathEdge;
        queue = other.queue;
        parentEdge = other.parentEdge;
        visitMark = other.visitMark;
        visitStamp = other.visitStamp;
        // The undo log is not copied: the field initialisers give the copy
        // an empty log with no trial running.
    }

    /** A copy that can change its capacities and flow independently. */
    public FlowNetwork copy() {
        return new FlowNetwork(this);
    }

    /**
     * Adds an edge u -> v with the given capacity (and its backward edge).
     *
     * @return the number of the forward edge
     */
    public int addEdge(int u, int v, int capacity) {
        int a = 2 * edgeCount++;
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

    /** The value of the current flow. */
    public int value() {
        return value;
    }

    /** The flow on the given edge. */
    public int flow(int edge) {
        return flow[edge];
    }

    /**
     * Changes the capacity of a forward edge. The caller must make sure the edge
     * does not carry more flow than the new capacity (see cancelOneUnit).
     */
    public void setCapacity(int edge, int capacity) {
        if (flow[edge] > capacity) {
            throw new IllegalStateException("edge carries more flow than the new capacity");
        }
        record(~edge, cap[edge]);
        cap[edge] = capacity;
    }

    /**
     * Runs Dinitz's algorithm from the current flow: while the sink can be
     * reached in the residual network, compute the layered residual network and add a
     * blocking flow (28 Oct, slide 231).
     *
     * @return the value of the maximum flow
     */
    public int maximise() {
        while (buildLayeredNetwork()) {
            System.arraycopy(head, 0, currentEdge, 0, nodeCount);
            int pushed;
            do {
                pushed = pushPathInLayeredNetwork();
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
    private boolean buildLayeredNetwork() {
        Arrays.fill(distance, -1);
        int front = 0;
        int back = 0;
        queue[back++] = source;
        distance[source] = 0;
        while (front < back) {
            int u = queue[front++];
            if (distance[sink] >= 0 && distance[u] >= distance[sink]) {
                break;
            }
            for (int a = head[u]; a >= 0; a = next[a]) {
                int v = to[a];
                if (distance[v] < 0 && cap[a] - flow[a] > 0) {
                    distance[v] = distance[u] + 1;
                    queue[back++] = v;
                }
            }
        }
        return distance[sink] >= 0;
    }

    /**
     * Finds one s-t path in the layered residual network and adds the
     * smallest residual capacity on it to the flow (slide 232). The search
     * is iterative (paths in a large grid can be long) and advances each
     * node's current edge past edges that are saturated or lead to a dead
     * end, so over one round every edge is skipped at most once. A dead end
     * is removed from the layered residual network.
     *
     * @return the amount pushed, or 0 if the flow is a blocking flow
     */
    private int pushPathInLayeredNetwork() {
        int depth = 0;
        int u = source;
        while (true) {
            if (u == sink) {
                int smallest = Integer.MAX_VALUE;
                for (int i = 0; i < depth; i++) {
                    int a = pathEdge[i];
                    smallest = Math.min(smallest, cap[a] - flow[a]);
                }
                for (int i = 0; i < depth; i++) {
                    push(pathEdge[i], smallest);
                }
                return smallest;
            }
            // Advance u's current edge to the next usable edge of the layered residual network.
            int a = currentEdge[u];
            while (a >= 0 && (distance[to[a]] != distance[u] + 1 || cap[a] - flow[a] <= 0)) {
                a = next[a];
            }
            currentEdge[u] = a;
            if (a >= 0) {
                pathEdge[depth++] = a;
                u = to[a];
            } else {
                // Dead end: remove u from the layered residual network and step back.
                if (u == source) {
                    return 0;
                }
                distance[u] = -1;
                int back = pathEdge[--depth];
                u = to[back ^ 1];
                currentEdge[u] = next[currentEdge[u]];
            }
        }
    }

    /**
     * The Ford/Fulkerson algorithm (22 Oct, slide 203) restricted to s-t
     * paths that start with the given edge out of the source.
     * This is enough after one source edge has been added to a network whose
     * flow was already maximum. Let X be the nodes the source could reach
     * before. No residual edge left X, otherwise the flow would not have been
     * maximum. An s-t path in the residual network that entered X could
     * therefore never leave it again, so every such path, now and after
     * each later round, starts with the new edge and stays outside X. Searching
     * from the head of that edge finds them all, and the breadth first search
     * stops as soon as it meets the sink, so it usually only looks at the
     * part of the grid around the new emitter instead of the whole network.
     *
     * @param startEdge a forward edge leaving the source
     * @return the value of the maximum flow
     */
    public int augmentFrom(int startEdge) {
        return augmentFrom(startEdge, null);
    }

    /**
     * As augmentFrom(startEdge), but the search may only visit nodes marked in
     * allowed (null allows all). With a restriction the result is a valid flow
     * that may be below the maximum; it is used for quick estimates.
     */
    public int augmentFrom(int startEdge, boolean[] allowed) {
        int start = to[startEdge];
        while (cap[startEdge] - flow[startEdge] > 0) {
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
                        parentEdge[v] = a;
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
            // Edges between grid nodes have capacity 1, so the smallest
            // residual capacity on the path is one unit.
            for (int v = sink; v != start; v = to[parentEdge[v] ^ 1]) {
                push(parentEdge[v], 1);
            }
            push(startEdge, 1);
            value += 1;
        }
        return value;
    }

    /**
     * The mirror image of augmentFrom: the Ford/Fulkerson algorithm restricted
     * to s-t paths that end with the given edge into the sink. The search
     * runs backwards from the tail of that edge to the source and stops when
     * it gets there.
     *
     * @param endEdge a forward edge entering the sink
     * @return true if at least one unit was pushed
     */
    public boolean augmentInto(int endEdge) {
        int end = to[endEdge ^ 1];
        boolean pushedAny = false;
        while (cap[endEdge] - flow[endEdge] > 0) {
            int stamp = ++visitStamp[0];
            int front = 0;
            int back = 0;
            queue[back++] = end;
            visitMark[end] = stamp;
            visitMark[sink] = stamp; // never walk back into the sink
            boolean found = false;
            while (front < back && !found) {
                int v = queue[front++];
                // Edge u -> v is usable if it has residual capacity; it is the
                // partner of the edge r stored at v.
                for (int r = head[v]; r >= 0; r = next[r]) {
                    int a = r ^ 1;
                    int u = to[r];
                    if (visitMark[u] != stamp && cap[a] - flow[a] > 0) {
                        visitMark[u] = stamp;
                        parentEdge[u] = a; // the edge leaving u towards the end
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
            int smallest = cap[endEdge] - flow[endEdge];
            for (int u = source; u != end; u = to[parentEdge[u]]) {
                int a = parentEdge[u];
                smallest = Math.min(smallest, cap[a] - flow[a]);
            }
            for (int u = source; u != end; u = to[parentEdge[u]]) {
                push(parentEdge[u], smallest);
            }
            push(endEdge, smallest);
            value += smallest;
            pushedAny = true;
        }
        return pushedAny;
    }

    /**
     * Removes one unit of flow that passes through the given forward edge,
     * keeping flow conservation everywhere else. The unit is traced forwards
     * to the sink and backwards to the source (or round a cycle back to the
     * edge) and every edge on the way loses one unit.
     * This is how the network reacts when a node that carries a pulse becomes
     * an emitter, or an emitter is removed: the affected pulse is cancelled
     * and the next maximise call is free to find a new route.
     *
     * @return the source edge the cancelled unit started from, or -1 if the
     *         unit went round a cycle
     */
    public int cancelOneUnit(int edge) {
        if (flow[edge] <= 0) {
            throw new IllegalStateException("edge carries no flow");
        }
        push(edge, -1);
        int tail = to[edge ^ 1];

        // Forwards: the head now receives one unit less than it sends.
        int v = to[edge];
        boolean closedCycle = false;
        while (v != sink) {
            if (v == tail) {
                closedCycle = true;
                break;
            }
            int out = outgoingFlowEdge(v);
            push(out, -1);
            v = to[out];
        }

        // Backwards: the tail now sends one unit less than it receives.
        if (closedCycle) {
            return -1;
        }
        int startEdge = edge;
        int u = tail;
        while (u != source) {
            startEdge = incomingFlowEdge(u);
            push(startEdge, -1);
            u = to[startEdge ^ 1];
        }
        value -= 1;
        return startEdge;
    }

    private void push(int edge, int amount) {
        record(edge, amount);
        flow[edge] += amount;
        flow[edge ^ 1] -= amount;
    }

    /**
     * Starts recording changes so that rollback can undo them. Trying out a
     * candidate emitter then costs only as much as the paths it changes,
     * instead of O(e) for copying the whole network. Trials cannot be nested.
     */
    public void beginTrial() {
        if (trailSize >= 0) {
            throw new IllegalStateException("a trial is already running");
        }
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
            int edge = trailEdge[i];
            if (edge >= 0) {
                flow[edge] -= trailValue[i];
                flow[edge ^ 1] += trailValue[i];
            } else {
                cap[~edge] = trailValue[i];
            }
        }
        trailSize = -1;
        value = valueBeforeTrial;
    }

    /**
     * Adds to mark every node at either end of an edge whose flow changed
     * since beginTrial, and returns the nodes that were newly marked so the
     * caller can clear them again.
     */
    public int[] markTrialNodes(boolean[] mark) {
        int[] added = new int[2 * Math.max(trailSize, 0)];
        int count = 0;
        for (int i = 0; i < trailSize; i++) {
            int edge = trailEdge[i];
            if (edge < 0) {
                continue; // a capacity change, not a flow change
            }
            for (int v : new int[] {to[edge], to[edge ^ 1]}) {
                if (!mark[v]) {
                    mark[v] = true;
                    added[count++] = v;
                }
            }
        }
        return Arrays.copyOf(added, count);
    }

    private void record(int edge, int amount) {
        if (trailSize < 0) {
            return; // not in a trial
        }
        if (trailSize == trailEdge.length) {
            trailEdge = Arrays.copyOf(trailEdge, 2 * trailSize);
            trailValue = Arrays.copyOf(trailValue, 2 * trailSize);
        }
        trailEdge[trailSize] = edge;
        trailValue[trailSize] = amount;
        trailSize++;
    }

    /** A forward edge leaving u that carries positive flow. */
    private int outgoingFlowEdge(int u) {
        for (int a = head[u]; a >= 0; a = next[a]) {
            if ((a & 1) == 0 && flow[a] > 0) {
                return a;
            }
        }
        throw new IllegalStateException("flow conservation violated");
    }

    /** A forward edge entering u that carries positive flow. */
    private int incomingFlowEdge(int u) {
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
     * (Max-Flow Min-Cut Theorem, 23 Oct, slides 215-219).
     */
    public boolean[] canReachSink() {
        return canReach(sink);
    }

    /**
     * Marks every node that can reach target in the residual network by a
     * path that does not pass through the sink (a backward breadth first
     * search from target).
     */
    public boolean[] canReach(int target) {
        boolean[] reach = new boolean[nodeCount];
        int front = 0;
        int back = 0;
        queue[back++] = target;
        reach[target] = true;
        while (front < back) {
            int v = queue[front++];
            // An edge u -> v is usable if it has residual capacity; its reverse
            // edge (stored at v) points back to u.
            for (int r = head[v]; r >= 0; r = next[r]) {
                int a = r ^ 1;
                int u = to[r];
                if (!reach[u] && u != sink && cap[a] - flow[a] > 0) {
                    reach[u] = true;
                    queue[back++] = u;
                }
            }
        }
        return reach;
    }
}

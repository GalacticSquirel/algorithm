package uk.ac.aber.cs31920.assignment;

/**
 * The energy grid: a square triangle mesh of side length l.
 *
 * Coordinates run from 0 to 2l in both directions and a pair (x, y) is a node
 * only if x and y are both even or both odd. That gives (l+1)^2 + l^2 nodes.
 *
 * Edges:
 * - an even node (x, y) is joined to the even nodes (x +- 2, y) and (x, y +- 2)
 *   along the square grid lines, and to the odd nodes (x +- 1, y +- 1) along the
 *   diagonals, so an interior even node has 8 neighbours, one on a side has 5
 *   and one in a corner has 3;
 * - an odd node (the centre of a square cell) is joined only to its four
 *   diagonal neighbours (x +- 1, y +- 1), so it always has 4 neighbours.
 *
 * Nodes are numbered 0 .. nodeCount()-1 so that the rest of the program can use
 * plain int arrays instead of coordinate objects.
 */
public final class Grid {

    /** Offsets to the even neighbours of an even node (along the grid lines). */
    private static final int[][] AXIS_STEPS = {{2, 0}, {-2, 0}, {0, 2}, {0, -2}};
    /** Offsets along the diagonals; used by both even and odd nodes. */
    private static final int[][] DIAGONAL_STEPS = {{1, 1}, {1, -1}, {-1, 1}, {-1, -1}};

    private final int sideLength;
    private final int maxCoordinate;
    /** idAt[x][y] is the node number of (x, y), or -1 if (x, y) is not a node. */
    private final int[][] idAt;
    private final int[] xOf;
    private final int[] yOf;
    /** neighbours[v] lists the node numbers adjacent to v. */
    private final int[][] neighbours;

    /**
     * Builds the grid of the given side length.
     *
     * @param sideLength the value l from the input file
     */
    public Grid(int sideLength) {
        this.sideLength = sideLength;
        this.maxCoordinate = 2 * sideLength;
        int nodes = (sideLength + 1) * (sideLength + 1) + sideLength * sideLength;
        idAt = new int[maxCoordinate + 1][maxCoordinate + 1];
        xOf = new int[nodes];
        yOf = new int[nodes];

        // Number the valid coordinate pairs row by row.
        int next = 0;
        for (int x = 0; x <= maxCoordinate; x++) {
            for (int y = 0; y <= maxCoordinate; y++) {
                if (x % 2 == y % 2) {
                    idAt[x][y] = next;
                    xOf[next] = x;
                    yOf[next] = y;
                    next++;
                } else {
                    idAt[x][y] = -1;
                }
            }
        }

        // Build the adjacency lists once; they never change.
        neighbours = new int[nodes][];
        int[] buffer = new int[8];
        for (int v = 0; v < nodes; v++) {
            int count = 0;
            if (xOf[v] % 2 == 0) {
                count = collect(v, AXIS_STEPS, buffer, count);
            }
            count = collect(v, DIAGONAL_STEPS, buffer, count);
            neighbours[v] = new int[count];
            System.arraycopy(buffer, 0, neighbours[v], 0, count);
        }
    }

    /** Appends to buffer every node reached from v by one of the given steps. */
    private int collect(int v, int[][] steps, int[] buffer, int count) {
        for (int[] step : steps) {
            int w = nodeAt(xOf[v] + step[0], yOf[v] + step[1]);
            if (w >= 0) {
                buffer[count++] = w;
            }
        }
        return count;
    }

    /** @return the node number of (x, y), or -1 if it is outside the grid or not a node */
    public int nodeAt(int x, int y) {
        if (x < 0 || y < 0 || x > maxCoordinate || y > maxCoordinate) {
            return -1;
        }
        return idAt[x][y];
    }

    public int nodeCount() {
        return xOf.length;
    }

    public int sideLength() {
        return sideLength;
    }

    public int x(int v) {
        return xOf[v];
    }

    public int y(int v) {
        return yOf[v];
    }

    public int[] neighbours(int v) {
        return neighbours[v];
    }

    /** The number of incident edges, i.e. the most pulses an emitter at v can send. */
    public int degree(int v) {
        return neighbours[v].length;
    }
}

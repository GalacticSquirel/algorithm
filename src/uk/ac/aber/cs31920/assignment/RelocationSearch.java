package uk.ac.aber.cs31920.assignment;

import java.util.ArrayList;
import java.util.List;

/**
 * Iterated local search (14 Oct, slides 108 and 110) for a working placement
 * with a fixed number k of emitters.
 *
 * A state is a set of k emitters, and its cost is the deficit: the pulses
 * still missing in the maximum flow. A move relocates one emitter, either to
 * a free node close to it or to one of the most promising free nodes for the
 * parasites that are still short. A move is kept if it does not make the
 * deficit worse; now and then a worse move is kept as well, so the search can
 * leave a plateau. After a run of moves without progress the search restarts
 * from the best state seen. Every state is checked with the exact maximum
 * flow, so a placement with deficit 0 really works.
 *
 * This is a heuristic: it can find a smaller placement quickly, but it can
 * never show that none exists. Proving that is left to the branch and bound
 * search in EmitterPlacer.
 */
public final class RelocationSearch {

    /** How far (in coordinate units) a "nearby" move may go. */
    private static final int NEARBY = 4;
    /** Moves without a new best state before restarting from the best state. */
    private static final int PATIENCE = 200;
    /** Promising candidates considered for a "targeted" move. */
    private static final int TARGETS = 8;

    private final Instance instance;
    private final Grid grid;
    private final List<Integer> forced;
    private long seed = 0x2545F4914F6CDD1DL;

    /**
     * @param instance the instance
     * @param forced   emitters that every solution contains; they never move
     */
    public RelocationSearch(Instance instance, List<Integer> forced) {
        this.instance = instance;
        this.grid = instance.grid();
        this.forced = forced;
    }

    /**
     * Looks for a working placement of k emitters, starting from start with
     * its least useful emitter dropped.
     *
     * @param start    a working placement with k + 1 emitters
     * @param maxMoves how many moves to try at most
     * @return a working placement with k emitters, or null if none was found
     */
    public List<Integer> search(List<Integer> start, int maxMoves) {
        PulseNetwork network = build(start);
        // Drop the emitter that currently sends the fewest pulses.
        int weakest = -1;
        for (int v : start) {
            if (!forced.contains(v) && (weakest < 0 || network.pulsesFrom(v) < network.pulsesFrom(weakest))) {
                weakest = v;
            }
        }
        if (weakest < 0) {
            return null; // only forced emitters, nothing can move
        }
        List<Integer> current = new ArrayList<>(start);
        current.remove(Integer.valueOf(weakest));
        network.removeEmitter(weakest);

        List<Integer> best = new ArrayList<>(current);
        int bestDeficit = network.deficit();
        int sinceBest = 0;
        for (int move = 0; move < maxMoves && network.deficit() > 0; move++) {
            if (sinceBest > PATIENCE) {
                current = new ArrayList<>(best);
                network = build(current);
                sinceBest = 0;
            }
            int index = pickMovable(current);
            if (index < 0) {
                return null;
            }
            int from = current.get(index);
            int to = (nextInt(3) == 0) ? nearby(from, network) : targeted(network);
            if (to < 0) {
                sinceBest++;
                continue;
            }
            PulseNetwork trial = network.copy();
            trial.removeEmitter(from);
            trial.addEmitter(to);
            if (trial.deficit() <= network.deficit() || nextInt(50) == 0) {
                network = trial;
                current.set(index, to);
            }
            if (network.deficit() < bestDeficit) {
                bestDeficit = network.deficit();
                best = new ArrayList<>(current);
                sinceBest = 0;
            } else {
                sinceBest++;
            }
        }
        return bestDeficit == 0 ? best : null;
    }

    /** A random emitter that is not forced, or -1 if there is none. */
    private int pickMovable(List<Integer> current) {
        int tries = 0;
        while (tries++ < 4 * current.size()) {
            int index = nextInt(current.size());
            if (!forced.contains(current.get(index))) {
                return index;
            }
        }
        return -1;
    }

    /** A random free node within NEARBY coordinate units of v, or -1. */
    private int nearby(int v, PulseNetwork network) {
        int x = grid.x(v) + nextInt(2 * NEARBY + 1) - NEARBY;
        int y = grid.y(v) + nextInt(2 * NEARBY + 1) - NEARBY;
        int w = grid.nodeAt(x, y);
        return (w >= 0 && network.isFree(w)) ? w : -1;
    }

    /**
     * One of the TARGETS free nodes with the highest estimated gain for the
     * current deficit, chosen at random, or -1 if no free node can help.
     */
    private int targeted(PulseNetwork network) {
        boolean[] reach = network.residualReachToSink();
        boolean[] useful = network.usefulNextEmitters(reach);
        int[] bestNodes = new int[TARGETS];
        int[] bestScores = new int[TARGETS];
        int count = 0;
        for (int v = 0; v < grid.nodeCount(); v++) {
            if (!useful[v]) {
                continue;
            }
            int score = network.promisingEdges(v, reach);
            // Keep the TARGETS highest scores (insertion into a small array).
            int pos;
            if (count < TARGETS) {
                pos = count++;
            } else if (score > bestScores[TARGETS - 1]) {
                pos = TARGETS - 1;
            } else {
                continue;
            }
            while (pos > 0 && bestScores[pos - 1] < score) {
                bestScores[pos] = bestScores[pos - 1];
                bestNodes[pos] = bestNodes[pos - 1];
                pos--;
            }
            bestScores[pos] = score;
            bestNodes[pos] = v;
        }
        return count == 0 ? -1 : bestNodes[nextInt(count)];
    }

    private PulseNetwork build(List<Integer> emitters) {
        PulseNetwork network = new PulseNetwork(instance);
        for (int v : emitters) {
            network.addEmitter(v);
        }
        return network;
    }

    /**
     * A pseudo-random number in 0 .. bound-1 from a xorshift generator. A
     * fixed seed keeps runs repeatable, so the program always prints the same
     * answer for the same input.
     */
    private int nextInt(int bound) {
        seed ^= seed << 13;
        seed ^= seed >>> 7;
        seed ^= seed << 17;
        return (int) Math.floorMod(seed, (long) bound);
    }
}

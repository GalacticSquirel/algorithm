package uk.ac.aber.cs31920.assignment;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

/**
 * A lower bound for crowded instances, based on regions that must contain an
 * emitter.
 *
 * Take a set R of free nodes and allow every free node outside R to be a
 * source (while still letting pulses pass through it). This relaxed network
 * allows at least every flow of every placement with no emitter in R. If its
 * maximum flow is below the total demand, then every working placement has an
 * emitter in R. Such an emitter sends at most cap(R) pulses, the largest
 * degree in R.
 *
 * Near a crowded corner or edge, R can be a few low-degree nodes, so the
 * emitter that has to be there wastes capacity. If R_1, ..., R_m are disjoint
 * regions that each must hold an emitter, then k emitters deliver at most
 * cap(R_1) + ... + cap(R_m) + 8 (k - m) pulses, which gives the bound
 *
 *     k >= m + ceil((total demand - cap(R_1) - ... - cap(R_m)) / 8).
 *
 * Forced emitters are regions of one node. The regions are found around each
 * parasite: start with the free nodes near it and, if they really must hold
 * an emitter, drop nodes (highest degree first) for as long as the rest still
 * must. Every test is a maximum flow (Dinitz, 28 Oct, slides 229-234, and
 * Ford-Fulkerson from one new source, 22 Oct, slide 203).
 */
public final class RegionBound {

    /** Largest radius (in coordinate units) of a region around a parasite. */
    private static final int MAX_RADIUS = 4;
    /** The radius whose region is also shrunk to its needed nodes. */
    private static final int SHRINK_RADIUS = 2;

    private final Instance instance;
    private final Grid grid;
    private final int maxEmitterPulses;

    public RegionBound(Instance instance, int maxEmitterPulses) {
        this.instance = instance;
        this.grid = instance.grid();
        this.maxEmitterPulses = maxEmitterPulses;
    }

    /** A set of free nodes that must hold at least count emitters. */
    private static final class Region {
        final int[] nodes;
        final int count;
        final int capacity;

        Region(int[] nodes, int count, int capacity) {
            this.nodes = nodes;
            this.count = count;
            this.capacity = capacity;
        }

        /** How much this region raises the bound over plain capacity counting. */
        int strength(int maxPulses) {
            return count * (maxPulses - capacity);
        }
    }

    /**
     * @param forced emitters every solution contains
     * @return a lower bound on the number of emitters
     */
    public int lowerBound(List<Integer> forced) {
        List<Region> regions = new ArrayList<>();
        for (int v : forced) {
            regions.add(new Region(new int[] {v}, 1, grid.degree(v)));
        }
        for (int p : instance.parasites()) {
            for (int radius = 1; radius <= MAX_RADIUS; radius++) {
                Region region = regionAround(p, radius, radius == SHRINK_RADIUS);
                if (region != null) {
                    regions.add(region);
                }
            }
        }
        // Pick disjoint regions greedily, strongest first.
        regions.sort(Comparator.comparingInt((Region r) -> -r.strength(maxEmitterPulses)));
        boolean[] used = new boolean[grid.nodeCount()];
        int count = 0;
        int capacity = 0;
        for (Region region : regions) {
            boolean disjoint = true;
            for (int v : region.nodes) {
                disjoint &= !used[v];
            }
            if (disjoint) {
                for (int v : region.nodes) {
                    used[v] = true;
                }
                count += region.count;
                capacity += region.count * region.capacity;
            }
        }
        int rest = Math.max(0, instance.totalDemand() - capacity);
        return count + (rest + maxEmitterPulses - 1) / maxEmitterPulses;
    }

    private int capacity(List<Integer> nodes) {
        int best = 0;
        for (int v : nodes) {
            best = Math.max(best, grid.degree(v));
        }
        return best;
    }

    /**
     * The free nodes within radius of p, if they must hold emitters: the
     * pulses that the relaxed network cannot bring in from outside have to
     * start inside, which needs ceil(missing / capacity) emitters there. With
     * shrink set, the region is then cut down (highest degree first) to the
     * nodes that are still needed, which lowers its capacity.
     *
     * @return the region, or null if it need not hold an emitter
     */
    private Region regionAround(int p, int radius, boolean shrink) {
        List<Integer> region = new ArrayList<>();
        boolean[] inRegion = new boolean[grid.nodeCount()];
        for (int dx = -radius; dx <= radius; dx++) {
            for (int dy = -radius; dy <= radius; dy++) {
                int v = grid.nodeAt(grid.x(p) + dx, grid.y(p) + dy);
                if (v >= 0 && !instance.hasParasite(v)) {
                    region.add(v);
                    inRegion[v] = true;
                }
            }
        }
        // Relaxed network: every free node outside the region is a source.
        PulseNetwork relaxed = new PulseNetwork(instance);
        for (int v = 0; v < grid.nodeCount(); v++) {
            if (!instance.hasParasite(v) && !inRegion[v]) {
                relaxed.openRelaxedSource(v);
            }
        }
        int missing = instance.totalDemand() - relaxed.maximiseRelaxed();
        if (missing <= 0) {
            return null;
        }
        if (!shrink) {
            int cap = capacity(region);
            return new Region(toArray(region), (missing + cap - 1) / cap, cap);
        }
        // Shrink: try opening the highest-degree nodes first.
        region.sort(Comparator.comparingInt((Integer v) -> -grid.degree(v)));
        List<Integer> kept = new ArrayList<>();
        for (int v : region) {
            relaxed.beginTrial();
            if (relaxed.augmentRelaxed(v) >= instance.totalDemand()) {
                relaxed.rollbackTrial(); // v is needed in the region
                kept.add(v);
            } else {
                relaxed.commitTrial(); // the rest still must hold an emitter
            }
        }
        return new Region(toArray(kept), 1, capacity(kept));
    }

    private static int[] toArray(List<Integer> list) {
        int[] result = new int[list.size()];
        for (int i = 0; i < result.length; i++) {
            result[i] = list.get(i);
        }
        return result;
    }
}

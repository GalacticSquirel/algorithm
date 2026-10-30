package uk.ac.aber.cs31920.assignment;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * One problem instance: the grid and the parasites on it.
 *
 * The input file holds the side length l on the first line and then one line
 * "x y demand" per parasite. The brief guarantees that files are well formed,
 * so the reader only skips blank lines and extra whitespace.
 */
public final class Instance {

    private final Grid grid;
    /** demand[v] is the demand of the parasite on node v, or 0 if v has no parasite. */
    private final int[] demand;
    /** The node numbers that hold a parasite, in input order. */
    private final List<Integer> parasites;
    private final int totalDemand;

    private Instance(Grid grid, int[] demand, List<Integer> parasites) {
        this.grid = grid;
        this.demand = demand;
        this.parasites = parasites;
        int sum = 0;
        for (int v : parasites) {
            sum += demand[v];
        }
        this.totalDemand = sum;
    }

    /**
     * Reads an instance from a file.
     *
     * @throws IOException if the file cannot be opened or read
     */
    public static Instance read(String fileName) throws IOException {
        List<int[]> numbers = new ArrayList<>();
        int sideLength = -1;
        try (BufferedReader in = new BufferedReader(new FileReader(fileName))) {
            String line;
            while ((line = in.readLine()) != null) {
                line = line.trim();
                if (line.isEmpty()) {
                    continue;
                }
                String[] parts = line.split("\\s+");
                if (sideLength < 0) {
                    sideLength = Integer.parseInt(parts[0]);
                } else {
                    numbers.add(new int[] {
                        Integer.parseInt(parts[0]), Integer.parseInt(parts[1]), Integer.parseInt(parts[2])
                    });
                }
            }
        }
        if (sideLength < 0) {
            throw new IOException("empty input file " + fileName);
        }

        Grid grid = new Grid(sideLength);
        int[] demand = new int[grid.nodeCount()];
        List<Integer> parasites = new ArrayList<>();
        for (int[] triple : numbers) {
            int v = grid.nodeAt(triple[0], triple[1]);
            if (v < 0) {
                throw new IOException("parasite at invalid node " + triple[0] + " " + triple[1]);
            }
            demand[v] = triple[2];
            parasites.add(v);
        }
        return new Instance(grid, demand, parasites);
    }

    public Grid grid() {
        return grid;
    }

    /** @return the demand of the parasite on v, or 0 if there is none */
    public int demand(int v) {
        return demand[v];
    }

    public boolean hasParasite(int v) {
        return demand[v] > 0;
    }

    public List<Integer> parasites() {
        return parasites;
    }

    public int totalDemand() {
        return totalDemand;
    }
}

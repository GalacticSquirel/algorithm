package uk.ac.aber.cs31920.assignment;

import java.io.IOException;
import java.util.List;

/**
 * Entry point. Every command line argument is an input file; for each one, in
 * order, a single line is printed: the emitter coordinates as "x y x y ...",
 * or "no solution" if the parasites cannot all be neutralised or the file
 * cannot be read.
 *
 * Usage: java uk.ac.aber.cs31920.assignment.MySolution input01.txt input02.txt
 */
public final class MySolution {

    /**
     * Time the improvement steps may spend on one input, in milliseconds. The
     * brief requires the fewest emitters, so there is no limit unless one is
     * asked for with -Dbudget=milliseconds (only meant for experiments).
     */
    private static final long SEARCH_BUDGET_MILLIS = Long.getLong("budget", Long.MAX_VALUE);

    private MySolution() {
    }

    public static void main(String[] args) {
        boolean verbose = System.getProperty("verbose") != null;
        for (String fileName : args) {
            long start = System.nanoTime();
            String line;
            String note = "";
            try {
                Instance instance = Instance.read(fileName);
                EmitterPlacer placer = new EmitterPlacer(instance, SEARCH_BUDGET_MILLIS);
                List<Integer> emitters = placer.solve();
                line = emitters == null ? "no solution" : format(instance.grid(), emitters);
                if (verbose && emitters != null) {
                    note = " [" + emitters.size() + " emitter(s), lower bound " + placer.lowerBound()
                        + (placer.provenOptimal() ? ", optimal" : ", not proven optimal") + "]";
                }
            } catch (IOException e) {
                // A file that cannot be opened or read counts as having no solution.
                line = "no solution";
                if (verbose) {
                    System.err.println(fileName + ": " + e.getMessage());
                }
            } catch (RuntimeException e) {
                // Should never happen. Report it on standard error, so that it
                // is not hidden, and still answer every remaining file.
                line = "no solution";
                e.printStackTrace();
            }
            System.out.println(line);
            if (verbose) {
                long millis = (System.nanoTime() - start) / 1_000_000;
                System.err.println(fileName + ": " + millis + " ms" + note);
            }
        }
    }

    /** Formats the emitters as "x y x y ...". */
    private static String format(Grid grid, List<Integer> emitters) {
        StringBuilder out = new StringBuilder();
        for (int v : emitters) {
            if (out.length() > 0) {
                out.append(' ');
            }
            out.append(grid.x(v)).append(' ').append(grid.y(v));
        }
        return out.toString();
    }
}

package uk.ac.aber.cs31920.assignment;

import java.io.IOException;
import java.util.List;

/**
 * Entry point. Every command line argument is an input file; for each one, in
 * order, a single line is printed: the emitter coordinates as "x y x y ...",
 * or "no solution" if the parasites cannot all be neutralised or the file
 * cannot be read.
 * Usage: java uk.ac.aber.cs31920.assignment.MySolution input01.txt input02.txt
 */
public final class MySolution {

    private MySolution() {
    }

    public static void main(String[] args) {
        for (String fileName : args) {
            String line;
            try {
                Instance instance = Instance.read(fileName);
                List<Integer> emitters = new EmitterPlacer(instance).solve();
                line = emitters == null ? "no solution" : format(instance.grid(), emitters);
            } catch (IOException e) {
                // A file that cannot be opened or read counts as having no solution.
                line = "no solution";
            } catch (RuntimeException e) {
                // Should never happen. Report it on standard error, so that it
                // is not hidden, and still answer every remaining file.
                line = "no solution";
                e.printStackTrace();
            }
            System.out.println(line);
        }
    }

    /** Formats the emitters as "x y x y ...". */
    private static String format(Grid grid, List<Integer> emitters) {
        StringBuilder out = new StringBuilder();
        for (int v : emitters) {
            if (!out.isEmpty()) {
                out.append(' ');
            }
            out.append(grid.x(v)).append(' ').append(grid.y(v));
        }
        return out.toString();
    }
}

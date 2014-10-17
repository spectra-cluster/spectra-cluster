package uk.ac.ebi.pride.spectracluster.util;

/**
 * CompareTo
 * User: Steve
 * Date: 2/5/14
 */
public class CompareTo {

    public static final double MINIMUM_DOUBLE_DIFFERENCE = 0.0000001;

    public static int compare(final double o1, final double o2) {
        double del = o1 - o2;
        double abs = Math.abs(del);
        if (abs < MINIMUM_DOUBLE_DIFFERENCE)
            return 0;
        return del < 0 ? -1 : 1;
    }


    public static int compare(final float o1, final float o2) {
        double del = o1 - o2;
        double abs = Math.abs(del);
        if (abs < MINIMUM_DOUBLE_DIFFERENCE)
            return 0;
        return del < 0 ? -1 : 1;
    }


    public static int compare(final long o1, final long o2) {
        if (o1 == o2)
            return 0;
        return o1 < o2 ? -1 : 1;
    }


    public static int compare(final int o1, final int o2) {
        if (o1 == o2)
            return 0;
        return o1 < o2 ? -1 : 1;
    }
}

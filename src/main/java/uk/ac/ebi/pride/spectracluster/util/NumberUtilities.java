package uk.ac.ebi.pride.spectracluster.util;

import java.text.NumberFormat;

/**
 * @author Rui Wang
 * @version $Id$
 */
public class NumberUtilities {
    private NumberUtilities() {
    }

    /**
     * convert a double into a String with a given precision
     * default double formatting is not very pretty
     *
     * @param in        int to convert
     * @param precision positive precision
     * @return non-null formatted string
     */
    public static String formatInt(int in, int precision) {
        StringBuilder ret = new StringBuilder(Integer.toString(in));
        if (ret.length() > precision)
            throw new IllegalArgumentException("Cannot write " + in + " in " + precision + " digits");
        while (ret.length() < precision)
            ret.insert(0, "0");
        return (ret.toString());
    }

    /**
     * convert a double into a String with a given precision
     * default double formatting is not very pretty
     *
     * @param in        non-null Double to convert
     * @param precision positive precision
     * @return non-null formatted string
     */
    public static String formatDouble(Double in, int precision) {
        return (formatDouble(in.doubleValue(), precision));
    }

    /**
     * convert a double into a String with a given precision
     * default double formatting is not very pretty
     *
     * @param in        double to convert
     * @param precision positive precision
     * @return non-null formatted string
     */
    public static String formatDouble(double in, int precision) {
        NumberFormat nf = NumberFormat.getInstance();
        nf.setMaximumFractionDigits(precision);
        return (nf.format(in));
    }

    /**
     * convert a double into a String with a given precision
     * default double formatting is not very pretty
     * This version defaults precision to 2
     *
     * @param in non-null Double to convert
     * @param in positive precision
     * @return non-null formatted string
     */
    public static String formatDouble(Double in) {
        return (formatDouble(in.doubleValue()));
    }

    /**
     * convert a double into a String with a given precision
     * default double formatting is not very pretty
     * This version defaults precision to 2
     *
     * @param in double to convert
     * @return non-null formatted string
     */
    public static String formatDouble(double in) {
        return (formatDouble(in, 2));
    }
}

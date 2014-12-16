package uk.ac.ebi.pride.spectracluster.filter;

import uk.ac.ebi.pride.spectracluster.spectrum.*;

import java.io.*;
import java.util.*;

/**
 * uk.ac.ebi.pride.spectracluster.filter.BinnedHighestNPeakFilter
 * return the highest peaks in binsize bins using overlapping bins -
 * default bin size is 100
 *
 * @author Steve Lewis
 * @date 27/05/2014
 */
@Deprecated // This code is scientifically highly problematic and should not be used
public class MaximialPeakFilter implements IPeakFilter {

    public static final int DEFAULT_MAX_PEAKS = 128;

    public static final int FIRST_BIN_MAX = 12; // todo This might be raised to filter less aggressively slewis

    /**
     * filters to be applied to cut peaks unter
     */
    public static final IPeakFilter[] DECREASING_BINS = {
            new BinnedHighestNPeakFilter(FIRST_BIN_MAX, BinnedHighestNPeakFilter.DEFAULT_BIN_SIZE, BinnedHighestNPeakFilter.DEFAULT_BIN_OVERLAP),
            new BinnedHighestNPeakFilter(FIRST_BIN_MAX - 1, BinnedHighestNPeakFilter.DEFAULT_BIN_SIZE, BinnedHighestNPeakFilter.DEFAULT_BIN_OVERLAP),
            new BinnedHighestNPeakFilter(FIRST_BIN_MAX - 2, BinnedHighestNPeakFilter.DEFAULT_BIN_SIZE, BinnedHighestNPeakFilter.DEFAULT_BIN_OVERLAP),
            new BinnedHighestNPeakFilter(FIRST_BIN_MAX - 3, BinnedHighestNPeakFilter.DEFAULT_BIN_SIZE, BinnedHighestNPeakFilter.DEFAULT_BIN_OVERLAP),
            new BinnedHighestNPeakFilter(FIRST_BIN_MAX - 4, BinnedHighestNPeakFilter.DEFAULT_BIN_SIZE, BinnedHighestNPeakFilter.DEFAULT_BIN_OVERLAP),
            new BinnedHighestNPeakFilter(FIRST_BIN_MAX - 5, BinnedHighestNPeakFilter.DEFAULT_BIN_SIZE, BinnedHighestNPeakFilter.DEFAULT_BIN_OVERLAP),
            new BinnedHighestNPeakFilter(FIRST_BIN_MAX - 6, BinnedHighestNPeakFilter.DEFAULT_BIN_SIZE, BinnedHighestNPeakFilter.DEFAULT_BIN_OVERLAP),
            new BinnedHighestNPeakFilter(FIRST_BIN_MAX - 7, BinnedHighestNPeakFilter.DEFAULT_BIN_SIZE, BinnedHighestNPeakFilter.DEFAULT_BIN_OVERLAP),
    };

    //============================================================
    /**
     * Statistice code - how are filters used
     */
    public static final int[] FILTER_USE_COUNTS = new int[DECREASING_BINS.length + 1];

    public static final int[] SPECTRUM_SIZE_COUNTS = new int[200];

    public static int NumberOverMax;

    public static void showStatistics(Appendable out) {
        try {
            int total = 0;
            final int[] sizeCounts = MaximialPeakFilter.SPECTRUM_SIZE_COUNTS;
            for (int i = 0; i < sizeCounts.length; i++) {
                int sizeCount = sizeCounts[i];
                total += sizeCount;
                out.append("" + 10 * i + " - " + 10 * (i + 1) + "\t" + sizeCount + "\n");
            }


            out.append("Total spectra " + total + "\n");
            out.append("Number unfiltered " + MaximialPeakFilter.NumberOverMax + "\n");
        } catch (IOException e) {
            throw new UnsupportedOperationException(e);
        }
    }

    //============================================================

    private final int maxPeaks;

    public MaximialPeakFilter(int maxPeaks) {
        this.maxPeaks = maxPeaks;
    }


    /**
     * Filter a given list of peaks
     *
     * @param peaks given list of peaks
     * @return a list of filtered peaks
     */
    @Override
    public List<IPeak> filter(List<IPeak> peaks) {
        List<IPeak> ret = peaks;
        int filterUsed = 0;

        int startSize = peaks.size();
        int startBin = Math.min(SPECTRUM_SIZE_COUNTS.length - 1, startSize / 10);

        SPECTRUM_SIZE_COUNTS[startBin]++; // count size distribution in bins of 10;

        // use incremental filters until the wanted number of peaks is reached
        while (ret.size() > maxPeaks) {
            ret = DECREASING_BINS[filterUsed++].filter(ret);
            if (filterUsed >= DECREASING_BINS.length) {
                NumberOverMax++; // we ran out of filters to apply
                break;
            }
        }

        FILTER_USE_COUNTS[filterUsed]++; // which filter did we use 0 is none
        return ret;
    }


}

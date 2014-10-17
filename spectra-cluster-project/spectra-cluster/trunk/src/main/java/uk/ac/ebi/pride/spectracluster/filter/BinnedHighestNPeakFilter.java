package uk.ac.ebi.pride.spectracluster.filter;

import uk.ac.ebi.pride.spectracluster.spectrum.IPeak;
import uk.ac.ebi.pride.spectracluster.util.MZIntensityUtilities;
import uk.ac.ebi.pride.spectracluster.util.comparator.PeakIntensityComparator;

import java.util.*;

/**
 * uk.ac.ebi.pride.spectracluster.filter.BinnedHighestNPeakFilter
 * return the highest peaks in binsize bins using overlapping bins -
 * default bin size is 100
 *
 * @author Steve Lewis
 * @date 27/05/2014
 */
public class BinnedHighestNPeakFilter implements IPeakFilter {

    public static final int MINIMUM_BINNED_MZ = MZIntensityUtilities.LOWEST_USABLE_MZ;
    public static final int MAXIMUM_BINNED_MZ = MZIntensityUtilities.HIGHEST_USABLE_MZ;
    public static final int DEFAULT_MAX_PEAKS_PER_BIN = 8;
    public static final int DEFAULT_BIN_SIZE = 100;
    public static final int DEFAULT_BIN_OVERLAP = DEFAULT_BIN_SIZE / 2;

    public static final IPeakFilter DEFAULT = new BinnedHighestNPeakFilter(DEFAULT_MAX_PEAKS_PER_BIN, DEFAULT_BIN_SIZE, DEFAULT_BIN_OVERLAP);

    public static final Comparator<IPeak> INTENSITY_COMPARATOR = PeakIntensityComparator.INSTANCE;

    private final int maxPeaks;
    private final int binSize;
    private final int binOverlap;

    public BinnedHighestNPeakFilter(int maxPeaks, int binSize, int binOverlap) {
        this.maxPeaks = maxPeaks;
        this.binSize = binSize;
        this.binOverlap = binOverlap;
    }

    public BinnedHighestNPeakFilter(int maxPeaks, int binSize) {
        this(maxPeaks, binSize, binSize / 2);
    }

    public BinnedHighestNPeakFilter(int maxPeaks) {
        this(maxPeaks, DEFAULT_BIN_SIZE);
    }

    public BinnedHighestNPeakFilter() {
        this(DEFAULT_MAX_PEAKS_PER_BIN);
    }


    /**
     * Filter a given list of peaks
     *
     * @param peaks given list of peaks
     * @return a list of filtered peaks
     */
    @Override
    public List<IPeak> filter(List<IPeak> peaks) {
        Set<IPeak> retained = new HashSet<IPeak>();
        int startpeak = 0;
        for (double binBottom = MINIMUM_BINNED_MZ; binBottom < MAXIMUM_BINNED_MZ - binSize; binBottom += binOverlap) {
            startpeak = handleBin(peaks, startpeak, retained, binBottom);
            if (startpeak > peaks.size())
                break;
        }
        List<IPeak> ret = new ArrayList<IPeak>(retained); // make a sorted list
        Collections.sort(ret); // sort by mz
        return ret;
    }

    /**
     * we generate
     *
     * @param allpeaks
     * @param startpeak
     * @param retained
     * @param binBottom
     * @return
     */
    protected int handleBin(List<IPeak> allpeaks, int startpeak, Set<IPeak> retained, double binBottom) {
        List<IPeak> byIntensity = new ArrayList<IPeak>();
        int nextBin = startpeak;
        double nextBinStartMZ = binBottom + binOverlap; // start of next bin
        double binEnd = binBottom + binSize; // end of this bin
        int index = startpeak;
        IPeak test = null;
        for (; index < allpeaks.size(); index++) {
            IPeak newTest = allpeaks.get(index);
            if(test != null)  {
                if(test.getMz() > newTest.getMz())  {   // out of order
                    if(Math.abs(test.getMz() - newTest.getMz()) > 1.2 * MZIntensityUtilities.SMALL_MZ_DIFFERENCE )
                        throw new IllegalStateException("Peaks are NOT Sorted by MZ");
                }
            }
            test = newTest;
            final float testMz = test.getMz();
            if (testMz < binBottom)
                continue;
            if (testMz < nextBinStartMZ)
                nextBin = index; // keep building next bin

            if (testMz > binEnd)
                break; // done with this bin
            byIntensity.add(test); // accumulate
        }

        // now sort highest intensity first
        Collections.sort(byIntensity, INTENSITY_COMPARATOR);
        // add highest maxPeaks to retained
        int numberAdded = 0;
        for (IPeak iPeak : byIntensity) {
            retained.add(iPeak);

            if (++numberAdded >= maxPeaks)
                break;
        }
        if (nextBin >= allpeaks.size())
            return Integer.MAX_VALUE; // finished all peaks - lets quit;
        return nextBin;

    }
}

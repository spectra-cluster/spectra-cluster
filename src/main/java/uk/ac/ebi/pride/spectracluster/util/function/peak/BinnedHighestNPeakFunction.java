package uk.ac.ebi.pride.spectracluster.util.function.peak;

import uk.ac.ebi.pride.spectracluster.spectrum.IPeak;
import uk.ac.ebi.pride.spectracluster.spectrum.Peak;
import uk.ac.ebi.pride.spectracluster.util.MZIntensityUtilities;
import uk.ac.ebi.pride.spectracluster.util.comparator.PeakIntensityComparator;
import uk.ac.ebi.pride.spectracluster.util.function.IFunction;

import java.util.*;

/**
 * Return the highest peaks in bin size bins using overlapping bins -
 * default bin size is 100
 *
 * @author Steve Lewis
 * @author Rui Wang
 * @version $Id$
 */
public class BinnedHighestNPeakFunction implements IFunction<List<IPeak>, List<IPeak>> {

    //public static final int MINIMUM_BINNED_MZ = MZIntensityUtilities.LOWEST_USABLE_MZ;
    public static final int MINIMUM_BINNED_MZ = 0; // this otherwise leads to quite unexpected behaviour
    public static final int MAXIMUM_BINNED_MZ = MZIntensityUtilities.HIGHEST_USABLE_MZ;
    public static final int DEFAULT_MAX_PEAKS_PER_BIN = 8;
    public static final int DEFAULT_BIN_SIZE = 100;

    public static final Comparator<IPeak> INTENSITY_COMPARATOR = PeakIntensityComparator.INSTANCE;

    /**
     * Number of peaks per bin
     */
    private final int maxPeaks;
    /**
     * Bin size in m/z
     */
    private final int binSize;
    /**
     * Overlap between two bins (default binSize / 2)
     */
    private final int binOverlap;

    /**
     * Creates a new BinnedHighestNPeakFilter
     * @param maxPeaks Maximum number of peaks per bin.
     * @param binSize Size of a bin in m/z
     * @param binOverlap Overlap between two bins in m/z
     */
    public BinnedHighestNPeakFunction(int maxPeaks, int binSize, int binOverlap) {
        this.maxPeaks = maxPeaks;
        this.binSize = binSize;
        this.binOverlap = binOverlap;

        if (binOverlap > binSize || binOverlap == binSize)
            throw new IllegalStateException("Bin overlap must be smaller than the bin size.");
    }

    public BinnedHighestNPeakFunction(int maxPeaks, int binSize) {
        this(maxPeaks, binSize, binSize / 2);
    }

    public BinnedHighestNPeakFunction(int maxPeaks) {
        this(maxPeaks, DEFAULT_BIN_SIZE);
    }

    public BinnedHighestNPeakFunction() {
        this(DEFAULT_MAX_PEAKS_PER_BIN);
    }

    @Override
    public List<IPeak> apply(List<IPeak> originalPeaks) {
        Set<IPeak> retained = new HashSet<IPeak>();
        int startPeak = 0;
        for (double binBottom = MINIMUM_BINNED_MZ; binBottom < MAXIMUM_BINNED_MZ - binSize; binBottom += (binSize - binOverlap)) {
            startPeak = handleBin(originalPeaks, startPeak, retained, binBottom);
            if (startPeak > originalPeaks.size())
                break;
        }

        // make a sorted list
        List<IPeak> ret = new ArrayList<IPeak>();
        for (IPeak iPeak : retained) {
            ret.add(new Peak(iPeak));
        }

        // sort by mz
        Collections.sort(ret);
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
        int startIndexNextBin = startpeak; // the index of the next bin's peak
        double binEnd = binBottom + binSize; // end of this bin
        double nextBinStartMZ = binEnd - binOverlap; // start of next bin

        // get all peaks within the current bin
        int index = startpeak;
        IPeak currentPeak = null;
        List<IPeak> byIntensity = new ArrayList<IPeak>();

        for (; index < allpeaks.size(); index++) {
            IPeak nextPeak = allpeaks.get(index);

            // make sure the peaks are sorted according to m/z
            if(currentPeak != null)  {
                if(currentPeak.getMz() > nextPeak.getMz())  {   // out of order
                    // only throw an exception if the difference is large enough
                    if(Math.abs(currentPeak.getMz() - nextPeak.getMz()) > 1.2 * MZIntensityUtilities.SMALL_MZ_DIFFERENCE )
                        throw new IllegalStateException("Peaks are NOT Sorted by MZ");
                }
            }

            // store all peaks that belong to this bin
            currentPeak = nextPeak;
            final float currentPeakMz = currentPeak.getMz();
            // ignore if it's before this bin
            if (currentPeakMz < binBottom)
                continue;
            // store all peaks that could belong to the next bin - iteratively building up to the right index
            if (currentPeakMz < nextBinStartMZ)
                startIndexNextBin = index;
            // done with this bin if we are above the end
            if (currentPeakMz > binEnd)
                break;

            byIntensity.add(currentPeak); // accumulate
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

        // finished all peaks - lets quit;
        if (startIndexNextBin >= allpeaks.size())
            return Integer.MAX_VALUE;


        return startIndexNextBin;
    }
}

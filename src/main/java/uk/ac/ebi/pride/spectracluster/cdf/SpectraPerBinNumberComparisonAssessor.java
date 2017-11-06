package uk.ac.ebi.pride.spectracluster.cdf;

import uk.ac.ebi.pride.spectracluster.cluster.ICluster;

/**
 * This INumberOfComparisonAssessor assess the number of
 * comparisons based on all spectra within a bin.
 *
 * Created by jg on 13.10.17.
 */
public class SpectraPerBinNumberComparisonAssessor implements INumberOfComparisonAssessor {
    private final float windowSize;
    private final int[] spectraPerBin;
    private final int maxBin;
    private final int minSpectra;
    public final static float MAX_PRECURSOR_MZ = 2500F;

    /**
     * Create a new SpectraPerBinNumberComparisonAssessor.
     * @param precursorTolerance The precursor tolerance in m/z to use.
     * @param minSpectra The minimum number of spectra always to return.
     */
    public SpectraPerBinNumberComparisonAssessor(float precursorTolerance, int minSpectra) {
        this.windowSize = precursorTolerance * 2;
        this.minSpectra = minSpectra;

        // initiate the bins
        int nBins = (int) Math.ceil(MAX_PRECURSOR_MZ / this.windowSize);
        maxBin = nBins - 1;
        spectraPerBin = new int[nBins];

        for (int i = 0; i < nBins; i++) {
            spectraPerBin[i] = 0;
        }
    }

    /**
     * Create a new SpectraPerBinNumberComparisonAssessor.
     * @param precursorTolerance The precursor tolerance in m/z to use.
     */
    public SpectraPerBinNumberComparisonAssessor(float precursorTolerance) {
        this(precursorTolerance, 0);
    }

    /**
     * Count this spectrum to know how many spectra exist per
     * bin. This function is thread safe.
     * @param precursorMz The spectrum's precursor m/z
     */
    public synchronized void countSpectrum(float precursorMz) {
        int bin = getBinForSpectrum(precursorMz);
        spectraPerBin[bin] += 1;
    }

    /**
     * Get the bin this spectrum belongs to.
     * @param precursorMz The spectrum's precursor m/z
     * @return 0-based bin.
     */
    private int getBinForSpectrum(float precursorMz) {
        int bin = (int) Math.floor(precursorMz / windowSize);
        if (bin > maxBin) {
            return maxBin;
        } else {
            return bin;
        }
    }

    @Override
    public int getNumberOfComparisons(ICluster clusterToCompare, int nCurrentClusters) {
        int bin = getBinForSpectrum(clusterToCompare.getPrecursorMz());

        int count;

        // never return anything lower than 1
        if (spectraPerBin[bin] <= 1)
            count = 1;
        else
            count = spectraPerBin[bin] - 1;

        if (count < minSpectra)
            return minSpectra;
        else
            return count;
    }
}

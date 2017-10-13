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
    public final static float MAX_PRECURSOR_MZ = 2500F;

    /**
     * Create a new SpectraPerBinNumberComparisonAssessor.
     * @param precursorTolerance The precursor tolerance in m/z to use.
     */
    public SpectraPerBinNumberComparisonAssessor(float precursorTolerance) {
        this.windowSize = precursorTolerance * 2;

        // initiate the bins
        int nBins = (int) Math.ceil(MAX_PRECURSOR_MZ / this.windowSize);
        maxBin = nBins - 1;
        spectraPerBin = new int[nBins];

        for (int i = 0; i < nBins; i++) {
            spectraPerBin[i] = 0;
        }
    }

    /**
     * Count this spectrum to know how many spectra exist per
     * bin.
     * @param precursorMz The spectrum's precursor m/z
     */
    public void countSpectrum(float precursorMz) {
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

        // never return anything lower than 1
        if (spectraPerBin[bin] <= 1)
            return 1;

        return spectraPerBin[bin] - 1;
    }
}

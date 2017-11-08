package uk.ac.ebi.pride.spectracluster.similarity;


import uk.ac.ebi.pride.spectracluster.spectrum.IPeak;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.util.Defaults;
import uk.ac.ebi.pride.spectracluster.util.Pair;


/**
 * Assesses the similarity between two
 * spectra using the normalized
 * dot-product as implemented by Frank et al. (2008) JPR
 * This implementation uses only the k highest peaks.
 * k is calculated by taking 15 peaks per 1000 Da
 * peptide mass. Furthermore, the vectors for the
 * dot-product are filled with the 1+ln(I) where I
 * is the peak's normalized intensity.
 * uk.ac.ebi.pride.spectracluster.similarity.FrankEtAlDotProduct
 *
 * @author jg
 */
public class FrankEtAlDotProduct implements ISimilarityChecker {
    private static final int K2011_BIN_SIZE = 50;

    /**
     * This is the minimum / default number of peaks to compare
     * as used in the original paper.
     */
    public final static int DEFAULT_NUMBER_OF_PEAKS_TO_COMPARE = 15;

    /**
     * If enabled the algorithm only uses the K
     * highest peaks of the spectra.
     */
    private boolean peakFiltering;

    /**
     * Peak filtering is disabled by default.
     */
    public final static boolean DEFAULT_PEAK_FILTERING = false;

    /**
     * The versions available from this algorithm. The only
     * difference is the way k is calculated.
     *
     * @author jg
     */
    public enum AlgorithmVersion {
        JPR_2008, NAT_METH_2011
    }

    public static final AlgorithmVersion DEFAULT_ALGORITHM = AlgorithmVersion.NAT_METH_2011;

    private float fragmentIonTolerance;
    private int numberOfPeaksToCompare;

    /**
     * Creates a new instance of the FrankEtAlDotProduct where all parameters
     * are set in the constructor.
     * @param fragmentIonTolerance The fragment ion tolerance to use for peak matching.
     * @param numberOfPeaksToCompare The default number of peaks to compare. This number is not used directly but
     *                               used within the two implementations (JPR_2008, NAT_METH_2011) of how the number
     *                               of peaks are being computed. The default value (15) is the one originally used
     *                               by the authors.
     * @param peakFiltering If enabled only the k highest peaks of each spectrum are being compared. If this is set
     *                      to false (default) the numberOfPeaksToCompare does not have any effect.
     */
    public FrankEtAlDotProduct(float fragmentIonTolerance,
                                  int numberOfPeaksToCompare,
                                  boolean peakFiltering) {
        this.fragmentIonTolerance = fragmentIonTolerance;
        this.numberOfPeaksToCompare = numberOfPeaksToCompare;
        this.peakFiltering = peakFiltering;
    }

    public FrankEtAlDotProduct(float fragmentIonTolerance) {
        this(fragmentIonTolerance, DEFAULT_NUMBER_OF_PEAKS_TO_COMPARE, DEFAULT_PEAK_FILTERING);
    }

    public FrankEtAlDotProduct() {
        this(Defaults.getFragmentIonTolerance(), DEFAULT_NUMBER_OF_PEAKS_TO_COMPARE, DEFAULT_PEAK_FILTERING);
    }

    /**
     * return a name which should not change
     *
     * @return !null name
     */
    @Override
    public String getName() {
        return getClass().getSimpleName();
    }

    /**
     * return a version number - this may be updated over time
     *
     * @return !null version
     */
    @Override
    public String getCurrentVersion() {
        return getVersion().toString();
    }

    /**
     * The algorithm version to use. By
     * default the version described in
     * Nature Methods 2011 will be used.
     */
    private AlgorithmVersion version = DEFAULT_ALGORITHM;

    @Override
    public double assessSimilarity(IPeakMatches peakMatches) {
        double dotProduct = 0;

        for (int i = 0; i < peakMatches.getNumberOfSharedPeaks(); i++) {
            Pair<IPeak, IPeak> matchedPeaks = peakMatches.getPeakPair(i);

            dotProduct += convertIntensity(matchedPeaks.getFirst()) * convertIntensity(matchedPeaks.getSecond());
        }

        // normalize the dot product
        double sumSquareIntensity1 = peakMatches.getSpectrumOne().getSumSquareIntensity();
        double sumSquareIntensity2 = peakMatches.getSpectrumTwo().getSumSquareIntensity();

        double denom = Math.sqrt(sumSquareIntensity1 * sumSquareIntensity2);
        if (denom == 0)
            return 0;
        double normalizedDotProduct = dotProduct / denom;

        if (normalizedDotProduct > 1.00000001) // JAVA rounding issue
            throw new IllegalStateException("Dot product must not exceed 1. (found " + normalizedDotProduct + ")");

        if (normalizedDotProduct > 1) // fix rounding issue
            normalizedDotProduct = 1;

        return normalizedDotProduct;
    }

    /**
     * Assesses the spectra's similarity using
     * the normalized dot-product
     */
    @Override
    public double assessSimilarity(ISpectrum spectrum1, ISpectrum spectrum2) {
        ISpectrum highestPeaksSpectrum1, highestPeaksSpectrum2;

        if (isPeakFiltering()) {
            // initialize the number of peaks1 to use with 15
            int numberCompared = computeNumberComparedSpectra(spectrum1, spectrum2);

            highestPeaksSpectrum1 = spectrum1.getHighestNPeaks(numberCompared);
            highestPeaksSpectrum2 = spectrum2.getHighestNPeaks(numberCompared);
        } else {
            // don't use peak filtering
            highestPeaksSpectrum1 = spectrum1;
            highestPeaksSpectrum2 = spectrum2;
        }

        IPeakMatches peakMatches = PeakMatchesUtilities.getSharedPeaksAsMatches(highestPeaksSpectrum1, highestPeaksSpectrum2, this.fragmentIonTolerance);

        return assessSimilarity(peakMatches);
    }

    /**
     * Transforms the intensities to penalize very high peaks.
     * This function is taken from the spectral-archives algorithm.
     */
    private double convertIntensity(IPeak p1) {
        double intensity = p1.getIntensity();
        if (intensity == 0)
            return 0;
        return 1 + Math.log(intensity);
    }

    protected int computeNumberComparedSpectra(ISpectrum spectrum1, ISpectrum spectrum2) {
        int numberComparedPeaks = numberOfPeaksToCompare;
        float precursorMz = spectrum1.getPrecursorMz();
        float precursor2 = spectrum2.getPrecursorMz();
        switch (version) {
            case NAT_METH_2011:
                numberComparedPeaks = calculateK2011(precursorMz, precursor2);
                break;
            case JPR_2008:
                numberComparedPeaks = calculateK2008(precursorMz, precursor2, spectrum1.getPrecursorCharge(), spectrum2.getPrecursorCharge());
                break;
        }
        return numberComparedPeaks;
    }


    /**
     * Calculate k by using 15 per 1000 Da of
     * peptide mass.
     *
     * @param precursor1
     * @param precursor2
     * @param charge1
     * @param charge2
     * @return
     */
    private int calculateK2008(Float precursor1, Float precursor2,
                               Integer charge1, Integer charge2) {
        // if any of the required values is missing, return 15
        if (precursor1 == null || precursor2 == null || charge1 == null || charge2 == null || charge1 <= 0 || charge2 <= 0)
            return numberOfPeaksToCompare;

        // take 15 peaks / 1000Da peptide mass
        double peptideMass = (precursor1 * charge1 + precursor2 * charge2) / 2;

        int largeBinningRegion = numberOfPeaksToCompare;
        int k = numberOfPeaksToCompare * (int) (peptideMass / largeBinningRegion);

        if (peptideMass % largeBinningRegion > 0)
            k += numberOfPeaksToCompare;

        return k;
    }

    /**
     * Calculate k by using the precursor m/z / 50.
     *
     * @param precursor1
     * @param precursor2
     * @return
     */
    private int calculateK2011(Float precursor1, Float precursor2) {
        // if any of the required values is missing, return 15
        if (precursor1 == null || precursor2 == null)
            return numberOfPeaksToCompare;

        // use m/z / 50

        return (int) ((precursor1 / K2011_BIN_SIZE + precursor2 / K2011_BIN_SIZE) / 2);
    }

    public AlgorithmVersion getVersion() {
        return version;
    }

    public void setVersion(AlgorithmVersion version) {
        this.version = version;
    }

    public boolean isPeakFiltering() {
        return peakFiltering;
    }

    public void setPeakFiltering(boolean peakFiltering) {
        this.peakFiltering = peakFiltering;
    }

    @Override
    public void setFragmentIonTolerance(float fragmentIonTolerance) {
        this.fragmentIonTolerance = fragmentIonTolerance;
    }

    @Override
    public float getFragmentIonTolerance() {
        return fragmentIonTolerance;
    }
}

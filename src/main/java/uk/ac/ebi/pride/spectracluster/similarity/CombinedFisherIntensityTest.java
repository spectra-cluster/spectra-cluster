package uk.ac.ebi.pride.spectracluster.similarity;

import org.apache.commons.math3.distribution.ChiSquaredDistribution;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.util.Defaults;

/**
 * This SimilarityChecker combines the probability of the
 * FisherExactTest and the IntensityRankCorrelation Test
 * using Fisher's method to combine extreme probabilities.
 *
 * Created by jg on 15.04.15.
 */
public class CombinedFisherIntensityTest implements ISimilarityChecker {
    public static final String algorithmName = "Combined FisherExact and Intensity rank test";
    public static final String algorithmVersion = "0.1";

    /**
     * These classes will be used to calculate the FisherExactTest and
     * IntensityRank probability.
     *
     * In both cases the IPeakMatches version of assessSimilarity
     * is used. Therefore, the fragmentIonTolerance and peakFiltering
     * do not have to be set.
      */
    protected final FisherExactTest fisherExactTest = new FisherExactTest();
    protected final IntensityRankCorrelation intensityRankCorrelation = new IntensityRankCorrelation();
    protected final ChiSquaredDistribution chiSquaredDistribution = new ChiSquaredDistribution(4); // always 4 degrees of freedom

    /**
     * The tolerance in m/z units used to match peaks
     */
    protected float fragmentIonTolerance;

    public static final boolean DEFAULT_PEAK_FILTERING = false;

    private boolean peakFiltering;

    public CombinedFisherIntensityTest() {
        this(Defaults.getFragmentIonTolerance());
    }

    public CombinedFisherIntensityTest(float fragmentIonTolerance) {
        this(fragmentIonTolerance, DEFAULT_PEAK_FILTERING);
    }

    public CombinedFisherIntensityTest(float fragmentIonTolerance, boolean peakFiltering) {
        this.fragmentIonTolerance = fragmentIonTolerance;
        this.peakFiltering = peakFiltering;
    }

    @Override
    public double assessSimilarity(ISpectrum spectrum1, ISpectrum spectrum2) {
        IPeakMatches peakMatches = PeakMatchesUtilities.getSharedPeaksAsMatches(spectrum1, spectrum2, fragmentIonTolerance, peakFiltering);
        return assessSimilarity(peakMatches);
    }

    @Override
    public double assessSimilarity(IPeakMatches peakMatches) {
        double fisherExactP = fisherExactTest.assessSimilarityAsPValue(peakMatches);
        double intensityRankP = intensityRankCorrelation.assessSimilarityAsPValue(peakMatches);

        // combine the p-values using Fisher's method
        double combined = -2 * (Math.log(fisherExactP) + Math.log(intensityRankP));
        double pValue;

        if (combined == 0)
            return 0;

        if (Double.isInfinite(combined))
            pValue = 0;
        else
            pValue = chiSquaredDistribution.density(combined); // TODO: this is not the correct implementation!

        return -Math.log(pValue);
    }

    @Override
    public boolean isPeakFiltering() {
        return peakFiltering;
    }

    @Override
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

    @Override
    public String getName() {
        return algorithmName;
    }

    @Override
    public String getCurrentVersion() {
        return algorithmVersion;
    }
}

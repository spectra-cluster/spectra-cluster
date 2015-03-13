package uk.ac.ebi.pride.spectracluster.similarity;

import org.apache.commons.math3.distribution.HypergeometricDistribution;

/**
 * This SimilarityChecker is based on the hypergeometric
 * probability that the observed similar m/z values are a
 * random match. It only calculates a point probability.
 *
 * Created by jg on 15.01.15.
 */
public class FisherExactTest extends HypergeometricScore {
    public static final String algorithmName = "Fisher Exact Test";
    public static final String algorithmVersion = "0.1";

    public FisherExactTest() {
        super();
    }

    public FisherExactTest(float peakMzTolerance) {
        super(peakMzTolerance);
    }

    public FisherExactTest(float peakMzTolerance, boolean peakFiltering) {
        super(peakMzTolerance, peakFiltering);
    }

    @Override
    protected double calculateSimilarityScore(int numberOfSharedPeaks, int numberOfPeaksFromSpec1, int numberOfPeaksFromSpec2, int numberOfBins) {
        HypergeometricDistribution hypergeometricDistribution = new HypergeometricDistribution(
                numberOfBins, numberOfPeaksFromSpec1, numberOfPeaksFromSpec2);

        // point probability
        double hgtScore = hypergeometricDistribution.probability(numberOfSharedPeaks);

        if (hgtScore == 0) {
            return 0;
        }

        return -Math.log(hgtScore);
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

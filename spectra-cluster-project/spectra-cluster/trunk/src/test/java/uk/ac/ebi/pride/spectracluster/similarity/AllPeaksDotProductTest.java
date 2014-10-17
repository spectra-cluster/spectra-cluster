package uk.ac.ebi.pride.spectracluster.similarity;

import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.util.ClusteringTestUtilities;
import uk.ac.ebi.pride.spectracluster.util.Defaults;

import java.util.List;

/**
 * Created by jg on 03.06.14.
 */
public class AllPeaksDotProductTest {
    private AllPeaksDotProduct allPeaksDotProduct;
    private List<ISpectrum> testSpectra;
    private FrankEtAlDotProduct frankEtAlDotProduct;

    @Before
    public void setUp() {
        allPeaksDotProduct = new AllPeaksDotProduct(Defaults.getSimilarityMZRange());
        frankEtAlDotProduct = new FrankEtAlDotProduct(Defaults.getSimilarityMZRange(), Defaults.getNumberComparedPeaks());
        testSpectra = ClusteringTestUtilities.readISpectraFromResource();
    }

    @Test
    public void testIdenticalSpectraComparison() {
        for (ISpectrum s : testSpectra) {
            double dotProduct = allPeaksDotProduct.assessSimilarity(s, s);
            Assert.assertEquals(1.0, dotProduct, 0);
        }
    }

    @Test
    public void testSimilarityChecker() {
        int totalComparisons = 0, acceptableRange = 0;

        // compare every spectrum with the next one
        for (int i = 0; i < testSpectra.size() - 1; i++) {
            ISpectrum current = testSpectra.get(i);
            ISpectrum next = testSpectra.get(i + 1);

            double dotProduct = allPeaksDotProduct.assessSimilarity(current, next);
            double originalDotProduct = frankEtAlDotProduct.assessSimilarity(current, next);

            totalComparisons++;

            if (Math.abs(originalDotProduct - dotProduct) <= 0.2) {
                acceptableRange++;
            }
        }

        double relativeAcceptable = (double) acceptableRange / (double) totalComparisons;
        System.out.println("total = " + totalComparisons + " acceptable = " + acceptableRange + " (" + relativeAcceptable + ")");
        Assert.assertTrue(relativeAcceptable >= 0.8);
    }
}

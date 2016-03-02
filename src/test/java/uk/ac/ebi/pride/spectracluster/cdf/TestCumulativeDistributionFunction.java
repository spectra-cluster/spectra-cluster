package uk.ac.ebi.pride.spectracluster.cdf;

import org.junit.Assert;
import org.junit.Test;
import uk.ac.ebi.pride.spectracluster.similarity.CombinedFisherIntensityTest;
import uk.ac.ebi.pride.spectracluster.similarity.FrankEtAlDotProduct;

/**
 * Created by jg on 05.05.15.
 */
public class TestCumulativeDistributionFunction {
    @Test
    public void testCombinedCdf() throws Exception {
        CumulativeDistributionFunction cumulativeDistributionFunction =
                CumulativeDistributionFunctionFactory.getCumulativeDistributionFunctionForSimilarityMetric(
                        CombinedFisherIntensityTest.class);

        Assert.assertEquals(6.965723908791688E-8, cumulativeDistributionFunction.probability(100, 4), 0.01);
        Assert.assertEquals(1.143539215064937E-5, cumulativeDistributionFunction.probability(60, 4), 0);

        Assert.assertTrue(cumulativeDistributionFunction.isSaveMatch(100, 4, 0.01));
        Assert.assertTrue(cumulativeDistributionFunction.isSaveMatch(60, 4, 0.01));
        Assert.assertFalse(cumulativeDistributionFunction.isSaveMatch(60, 40000000, 0.01));
        Assert.assertFalse(cumulativeDistributionFunction.isSaveMatch(30, 4000, 0.01));
    }

    @Test
    public void testDotCDf() throws Exception {
        CumulativeDistributionFunction cdf = CumulativeDistributionFunctionFactory.getCumulativeDistributionFunctionForSimilarityMetric(FrankEtAlDotProduct.class);

        Assert.assertEquals(1.5266607675812338E-6, cdf.probability(0.7, 4), 0.01);
        Assert.assertEquals(0.07722590845447774, cdf.probability(0.2, 10), 0.001);

        Assert.assertEquals(0.915328240919341, cdf.probability(0, 10), 0.001);
        Assert.assertEquals(0.0, cdf.probability(1, 10), 0.001);
        Assert.assertEquals(0.0, cdf.probability(1.1, 10), 0.001);
        Assert.assertEquals(0.0011151285305012193, cdf.probability(0.44, 10), 0.001);
    }
}

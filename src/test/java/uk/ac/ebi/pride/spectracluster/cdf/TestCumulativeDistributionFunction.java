package uk.ac.ebi.pride.spectracluster.cdf;

import org.junit.Assert;
import org.junit.Test;
import uk.ac.ebi.pride.spectracluster.similarity.CombinedFisherIntensityTest;

/**
 * Created by jg on 05.05.15.
 */
public class TestCumulativeDistributionFunction {
    @Test
    public void testCombinedCdf() throws Exception {
        CumulativeDistributionFunction cumulativeDistributionFunction =
                CumulativeDistributionFunctionFactory.getCumulativeDistributionFunctionForSimilarityMetric(
                        CombinedFisherIntensityTest.class);

        Assert.assertEquals(0.0, cumulativeDistributionFunction.probability(100, 4), 0);
        Assert.assertEquals(8.385432087587574E-7, cumulativeDistributionFunction.probability(60, 4), 0);

        Assert.assertTrue(cumulativeDistributionFunction.isSaveMatch(100, 4, 0.01));
        Assert.assertTrue(cumulativeDistributionFunction.isSaveMatch(60, 4, 0.01));
        Assert.assertFalse(cumulativeDistributionFunction.isSaveMatch(60, 40000000, 0.01));
        Assert.assertFalse(cumulativeDistributionFunction.isSaveMatch(30, 4000, 0.01));
    }
}

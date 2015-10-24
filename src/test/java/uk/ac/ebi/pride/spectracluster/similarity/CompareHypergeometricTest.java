package uk.ac.ebi.pride.spectracluster.similarity;

import cern.jet.random.HyperGeometric;
import cern.jet.random.Normal;
import cern.jet.random.engine.RandomEngine;
import org.apache.commons.math3.distribution.HypergeometricDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.junit.Assert;
import org.junit.Test;

/**
 * Created by jg on 21.05.15.
 */
public class CompareHypergeometricTest {
    @Test
    public void compareHypergeometric() throws Exception {
        int nBins = 200;
        int binsFromSpec1 = 30;
        int binsFromSpec2 = 40;
        int sharedPeaks = 5;

        double d1 = new HypergeometricDistribution(nBins, binsFromSpec1, binsFromSpec2).probability(sharedPeaks);
        double d2 = new HyperGeometric(nBins, binsFromSpec1, binsFromSpec2, RandomEngine.makeDefault()).pdf(sharedPeaks);

        Assert.assertEquals(d1, d2, 0.00001);
    }

    @Test
    public void compareNormal() throws Exception {
        double sd = 1;
        double correlation = 0.4;


        NormalDistribution correlationDistribution = new NormalDistribution(0, sd);
        double probability1 = correlationDistribution.cumulativeProbability(correlation);

        Normal normal = new Normal(0, sd, FisherExactTest.randomEngine);
        double probability2 = normal.cdf(correlation);

        Assert.assertEquals(probability1, probability2, 0.000001);
    }
}

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
        allPeaksDotProduct = new AllPeaksDotProduct(Defaults.getFragmentIonTolerance());
        frankEtAlDotProduct = new FrankEtAlDotProduct(Defaults.getFragmentIonTolerance(), Defaults.getNumberComparedPeaks());
        testSpectra = ClusteringTestUtilities.readISpectraFromResource();
    }

    @Test
    public void testIdenticalSpectraComparison() {
        for (ISpectrum s : testSpectra) {
            double dotProduct = allPeaksDotProduct.assessSimilarity(s, s);
            Assert.assertEquals(1.0, dotProduct, 0);
        }
    }
}

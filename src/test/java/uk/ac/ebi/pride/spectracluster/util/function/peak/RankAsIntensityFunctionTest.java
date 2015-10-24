package uk.ac.ebi.pride.spectracluster.util.function.peak;

import junit.framework.Assert;
import org.junit.Test;
import uk.ac.ebi.pride.spectracluster.spectrum.IPeak;
import uk.ac.ebi.pride.spectracluster.spectrum.Peak;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by jg on 18.04.15.
 */
public class RankAsIntensityFunctionTest {
    private RankAsIntensityFunction rankAsIntensityFunction = new RankAsIntensityFunction();

    @Test
    public void testIntensityConversion() {
        List<IPeak> testPeakList = new ArrayList<IPeak>(5);

        IPeak p1 = new Peak(1F, 20F);
        IPeak p2 = new Peak(2F, 30F);
        IPeak p3 = new Peak(3F, 40F);
        IPeak p4 = new Peak(4F, 50F);

        testPeakList.add(p3);
        testPeakList.add(p2);
        testPeakList.add(p1);
        testPeakList.add(p4);

        List<IPeak> convertedPeaks = rankAsIntensityFunction.apply(testPeakList);

        Assert.assertEquals(3F, convertedPeaks.get(0).getIntensity());
        Assert.assertEquals(2F, convertedPeaks.get(1).getIntensity());
        Assert.assertEquals(1F, convertedPeaks.get(2).getIntensity());
        Assert.assertEquals(4F, convertedPeaks.get(3).getIntensity());
    }
}

package uk.ac.ebi.pride.spectracluster.util.function.spectrum;

import junit.framework.Assert;
import org.junit.Before;
import org.junit.Test;
import uk.ac.ebi.pride.spectracluster.spectrum.IPeak;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.spectrum.Peak;
import uk.ac.ebi.pride.spectracluster.spectrum.Spectrum;
import uk.ac.ebi.pride.spectracluster.util.Defaults;
import uk.ac.ebi.pride.spectracluster.util.function.IFunction;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by jg on 06.11.17.
 */
public class BinSpectrumFunctionTest {
    private ISpectrum testSpectrum;

    @Before
    public void setUp() {
        List<IPeak> peakList = new ArrayList<IPeak>(5);

        peakList.add(new Peak(100, 20, 1));
        peakList.add(new Peak(100.2F, 5, 1));
        peakList.add(new Peak(100.1F, 10, 1));
        peakList.add(new Peak(103, 20, 1));
        peakList.add(new Peak(103.4F, 10, 3));

        testSpectrum = new Spectrum("test", 2, 10, Defaults.getDefaultQualityScorer(), peakList);
    }

    @Test
    public void testAvPeakBinning() throws Exception {
        IFunction<ISpectrum, ISpectrum> binFunction = new BinSpectrumAverageFunction(0.5F);

        ISpectrum binnedSpectrum = binFunction.apply(testSpectrum);
        List<IPeak> binnedPeaks = binnedSpectrum.getPeaks();

        Assert.assertEquals(2, binnedPeaks.size());
        Assert.assertEquals(35.0F, binnedPeaks.get(0).getIntensity());
        Assert.assertEquals(1, binnedPeaks.get(0).getCount());
        Assert.assertEquals(100.05F, binnedPeaks.get(0).getMz(), 0.01F);

        Assert.assertEquals(30.0F, binnedPeaks.get(1).getIntensity());
        Assert.assertEquals(3, binnedPeaks.get(1).getCount());
        Assert.assertEquals(103.13F, binnedPeaks.get(1).getMz(), 0.01F);
    }

    @Test
    public void testMaxPeakBinning() throws Exception {
        IFunction<ISpectrum, ISpectrum> binFunction = new BinSpectrumMaxFunction(0.5F);

        ISpectrum binnedSpectrum = binFunction.apply(testSpectrum);
        List<IPeak> binnedPeaks = binnedSpectrum.getPeaks();

        Assert.assertEquals(2, binnedPeaks.size());
        Assert.assertEquals(20.0F, binnedPeaks.get(0).getIntensity());
        Assert.assertEquals(1, binnedPeaks.get(0).getCount());
        Assert.assertEquals(100.0F, binnedPeaks.get(0).getMz(), 0.01F);

        Assert.assertEquals(20.0F, binnedPeaks.get(1).getIntensity());
        Assert.assertEquals(1, binnedPeaks.get(1).getCount());
        Assert.assertEquals(103.0F, binnedPeaks.get(1).getMz(), 0.01F);
    }
}

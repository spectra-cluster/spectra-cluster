package uk.ac.ebi.pride.spectracluster.util.function.spectrum;

import junit.framework.Assert;
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
 * Created by jg on 05.07.16.
 */
public class RemovePrecursorPeaksFunctionTest {
    @Test
    public void testPrecursorRemoval() throws Exception {
        List<IPeak> peakList = new ArrayList<IPeak>();
        // 800.4702 808.9828 809.4834
        peakList.add(new Peak(800.4702f, 1));
        peakList.add(new Peak(808.9828f, 1));
        peakList.add(new Peak(809.4834f, 1));
        ISpectrum spectrum = new Spectrum("id1", 2, 808.9858f, Defaults.getDefaultQualityScorer(), peakList);

        IFunction<ISpectrum, ISpectrum> filter = new RemovePrecursorPeaksFunction(0.05f);

        ISpectrum filteredSpectrum = filter.apply(spectrum);

        Assert.assertEquals(0, filteredSpectrum.getPeaks().size());
    }
}

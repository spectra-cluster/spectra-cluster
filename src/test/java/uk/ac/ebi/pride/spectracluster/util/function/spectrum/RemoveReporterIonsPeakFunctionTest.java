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
public class RemoveReporterIonsPeakFunctionTest {
    @Test
    public void testRemoveReporterIons() throws Exception {
        List<IPeak> peakList = new ArrayList<IPeak>();
        // 130.1348 130.1410 131.1380
        peakList.add(new Peak(130.1348f, 1));
        peakList.add(new Peak(130.1410f, 1));
        peakList.add(new Peak(131.1380f, 1));
        ISpectrum spectrum = new Spectrum("id1", 2, 808.9858f, Defaults.getDefaultQualityScorer(), peakList);

        IFunction<ISpectrum, ISpectrum> filter = new RemoveReporterIonPeaksFunction(0.05f, RemoveReporterIonPeaksFunction.REPORTER_TYPE.TMT);

        ISpectrum filteredSpectrum = filter.apply(spectrum);

        Assert.assertEquals(0, filteredSpectrum.getPeaks().size());
    }
}

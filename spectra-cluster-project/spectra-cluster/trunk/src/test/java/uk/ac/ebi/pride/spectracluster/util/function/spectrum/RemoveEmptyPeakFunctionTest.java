package uk.ac.ebi.pride.spectracluster.util.function.spectrum;

import org.junit.Test;
import uk.ac.ebi.pride.spectracluster.spectrum.IPeak;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.spectrum.Peak;

import java.util.ArrayList;
import java.util.List;

import static org.junit.Assert.assertEquals;
import static org.mockito.Mockito.mock;
import static org.mockito.Mockito.when;

public class RemoveEmptyPeakFunctionTest {

    @Test
    public void testApply() throws Exception {
        RemoveEmptyPeakFunction removeEmptyPeakFunction = new RemoveEmptyPeakFunction();

        ISpectrum spectrum = mock(ISpectrum.class);
        List<IPeak> peaks = new ArrayList<IPeak>();
        peaks.add(new Peak(10f, 10f));
        peaks.add(new Peak(0f, 10f));
        peaks.add(new Peak(10f, 0f));
        peaks.add(new Peak(0f, 0f));
        when(spectrum.getPeaks()).thenReturn(peaks);

        ISpectrum newSpectrum = removeEmptyPeakFunction.apply(spectrum);
        assertEquals(1, newSpectrum.getPeaks().size());
    }
}
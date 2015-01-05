package uk.ac.ebi.pride.spectracluster.util.predicate.spectrum;

import org.junit.Test;
import uk.ac.ebi.pride.spectracluster.spectrum.IPeak;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.spectrum.Peak;

import java.util.ArrayList;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import static org.mockito.Mockito.mock;
import static org.mockito.Mockito.when;

public class NoneNegativePeakPredicateTest {

    @Test
    public void testApply() throws Exception {
        NoneNegativePeakPredicate noneNegativePeakPredicate = new NoneNegativePeakPredicate();

        ISpectrum spectrum = mock(ISpectrum.class);
        ArrayList<IPeak> peaks = new ArrayList<IPeak>();

        peaks.add(new Peak(-100, 100));
        when(spectrum.getPeaks()).thenReturn(peaks);

        assertFalse(noneNegativePeakPredicate.apply(spectrum));

        peaks.clear();
        peaks.add(new Peak(100, -100));
        assertFalse(noneNegativePeakPredicate.apply(spectrum));

        peaks.clear();
        peaks.add(new Peak(100, 100));
        assertTrue(noneNegativePeakPredicate.apply(spectrum));
    }
}
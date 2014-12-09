package uk.ac.ebi.pride.spectracluster.util.predicate.spectrum;

import org.junit.Test;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import static org.mockito.Mockito.mock;
import static org.mockito.Mockito.when;

public class MinimumNumberOfPeaksPredicateTest {

    @Test
    public void testApply() throws Exception {
        MinimumNumberOfPeaksPredicate MinimumNumberOfPeaksPredicate = new MinimumNumberOfPeaksPredicate();

        ISpectrum spectrum = mock(ISpectrum.class);
        when(spectrum.getPeaksCount()).thenReturn(200);
        assertTrue(MinimumNumberOfPeaksPredicate.apply(spectrum));

        when(spectrum.getPeaksCount()).thenReturn(100);
        assertTrue(MinimumNumberOfPeaksPredicate.apply(spectrum));

        when(spectrum.getPeaksCount()).thenReturn(10);
        assertFalse(MinimumNumberOfPeaksPredicate.apply(spectrum));
    }
}
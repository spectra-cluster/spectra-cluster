package uk.ac.ebi.pride.spectracluster.util.predicate.spectrum;

import org.junit.Test;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import static org.mockito.Mockito.mock;
import static org.mockito.Mockito.when;

public class WithPrecursorPredicateTest {

    @Test
    public void testApply() throws Exception {
        WithPrecursorPredicate withPrecursorPredicate = new WithPrecursorPredicate();

        ISpectrum spectrum = mock(ISpectrum.class);
        when(spectrum.getPrecursorMz()).thenReturn(200f);
        when(spectrum.getPrecursorCharge()).thenReturn(-1);
        assertTrue(withPrecursorPredicate.apply(spectrum));

        when(spectrum.getPrecursorMz()).thenReturn(1f);
        assertFalse(withPrecursorPredicate.apply(spectrum));

        when(spectrum.getPrecursorMz()).thenReturn(200f);
        when(spectrum.getPrecursorCharge()).thenReturn(0);
        assertFalse(withPrecursorPredicate.apply(spectrum));
    }
}
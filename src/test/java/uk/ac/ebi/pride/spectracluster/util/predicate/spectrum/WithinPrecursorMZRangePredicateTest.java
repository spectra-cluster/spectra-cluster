package uk.ac.ebi.pride.spectracluster.util.predicate.spectrum;

import org.junit.Test;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import static org.mockito.Mockito.mock;
import static org.mockito.Mockito.when;

public class WithinPrecursorMZRangePredicateTest {

    @Test
    public void testApply() throws Exception {
        WithinPrecursorMZRangePredicate withinPrecursorMZRangePredicate = new WithinPrecursorMZRangePredicate(10f, 200f);

        ISpectrum spectrum = mock(ISpectrum.class);
        when(spectrum.getPrecursorMz()).thenReturn(200f);
        assertTrue(withinPrecursorMZRangePredicate.apply(spectrum));

        when(spectrum.getPrecursorMz()).thenReturn(10f);
        assertTrue(withinPrecursorMZRangePredicate.apply(spectrum));

        when(spectrum.getPrecursorMz()).thenReturn(9f);
        assertFalse(withinPrecursorMZRangePredicate.apply(spectrum));

        when(spectrum.getPrecursorMz()).thenReturn(210f);
        assertFalse(withinPrecursorMZRangePredicate.apply(spectrum));

    }
}
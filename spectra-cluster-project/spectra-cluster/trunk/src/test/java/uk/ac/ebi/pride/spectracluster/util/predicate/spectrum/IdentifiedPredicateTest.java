package uk.ac.ebi.pride.spectracluster.util.predicate.spectrum;

import org.junit.Test;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.spectrum.KnownProperties;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import static org.mockito.Mockito.mock;
import static org.mockito.Mockito.when;

public class IdentifiedPredicateTest {

    @Test
    public void testApply() throws Exception {
        IdentifiedPredicate identifiedPredicate = new IdentifiedPredicate();

        ISpectrum spectrum = mock(ISpectrum.class);
        when(spectrum.getProperty(KnownProperties.IDENTIFIED_PEPTIDE_KEY)).thenReturn("KKKKKK");
        assertTrue(identifiedPredicate.apply(spectrum));

        when(spectrum.getProperty(KnownProperties.IDENTIFIED_PEPTIDE_KEY)).thenReturn(null);
        assertFalse(identifiedPredicate.apply(spectrum));
    }
}
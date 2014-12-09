package uk.ac.ebi.pride.spectracluster.util.predicate.spectrum;

import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.util.predicate.IPredicate;

/**
 * Check whether spectrum precursor details are present
 *
 * @author Rui Wang
 * @version $Id$
 */
public class WithPrecursorPredicate implements IPredicate<ISpectrum> {
    public static final float MINIMUM_PRECURSOR_MZ = 10f;

    @Override
    public boolean apply(ISpectrum spectrum) {
        return spectrum.getPrecursorMz() >= MINIMUM_PRECURSOR_MZ && spectrum.getPrecursorCharge() != 0;
    }
}

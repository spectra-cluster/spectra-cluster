package uk.ac.ebi.pride.spectracluster.util.predicate.spectrum;

import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.spectrum.KnownProperties;
import uk.ac.ebi.pride.spectracluster.util.predicate.IPredicate;

/**
 * Check whether a spectrum is unidentified or not
 *
 * @author Rui Wang
 * @version $Id$
 */
public class UnidentifiedPredicate implements IPredicate<ISpectrum> {

    @Override
    public boolean apply(ISpectrum spectrum) {
        return spectrum.getProperty(KnownProperties.IDENTIFIED_PEPTIDE_KEY) == null;
    }
}
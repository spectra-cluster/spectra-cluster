package uk.ac.ebi.pride.spectracluster.util.predicate.spectrum;

import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.spectrum.KnownProperties;
import uk.ac.ebi.pride.spectracluster.util.predicate.IPredicate;

/**
 * Check whether a spectrum is from a given taxonomy or not
 *
 * @author Rui Wang
 * @version $Id$
 */
public class TaxonomyPredicate implements IPredicate<ISpectrum> {

    private final String taxonomy;

    public TaxonomyPredicate(String taxonomy) {
        this.taxonomy = taxonomy;
    }

    @Override
    public boolean apply(ISpectrum spectrum) {
        String tax = spectrum.getProperty(KnownProperties.TAXONOMY_KEY);
        return tax != null && tax.equals(taxonomy);
    }


}

package uk.ac.ebi.pride.spectracluster.util.predicate.spectrum;

import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.util.predicate.IPredicate;

/**
 * Detect whether a spectrum's precursor m/z is within a given range
 *
 * @author Rui Wang
 * @version $Id$
 */
public class WithinPrecursorMZRangePredicate implements IPredicate<ISpectrum> {

    private final float miniMZ;
    private final float maxMZ;

    public WithinPrecursorMZRangePredicate(float miniMZ, float maxMZ) {
        if(miniMZ < 0 || miniMZ >= maxMZ)
            throw new IllegalArgumentException("bad MZ Limits min " + miniMZ + " max " + maxMZ);

        this.miniMZ = miniMZ;
        this.maxMZ = maxMZ;
    }

    @Override
    public boolean apply(ISpectrum spectrum) {
        float precursorMz = spectrum.getPrecursorMz();

        return precursorMz >= miniMZ && precursorMz <= maxMZ;
    }
}

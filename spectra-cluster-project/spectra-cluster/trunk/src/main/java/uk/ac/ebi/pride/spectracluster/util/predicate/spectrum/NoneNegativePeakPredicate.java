package uk.ac.ebi.pride.spectracluster.util.predicate.spectrum;

import uk.ac.ebi.pride.spectracluster.spectrum.IPeak;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.util.predicate.IPredicate;

import java.util.List;

/**
 * Check whether a given spectrum contains peaks that have negative m/z or intensity values
 *
 * @author Rui Wang
 * @version $Id$
 */
public class NoneNegativePeakPredicate implements IPredicate<ISpectrum> {
    @Override
    public boolean apply(ISpectrum spectrum) {
        List<IPeak> peaks = spectrum.getPeaks();

        for (IPeak peak : peaks) {
            if (peak.getMz() < 0 || peak.getIntensity() < 0) {
                return false;
            }
        }
        return true;
    }
}

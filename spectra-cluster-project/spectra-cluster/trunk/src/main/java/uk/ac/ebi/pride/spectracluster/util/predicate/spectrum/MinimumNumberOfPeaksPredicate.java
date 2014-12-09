package uk.ac.ebi.pride.spectracluster.util.predicate.spectrum;

import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.util.predicate.IPredicate;

/**
 * Check whether a spectrum has minimum number of peaks
 *
 * @author Rui Wang
 * @version $Id$
 */
public class MinimumNumberOfPeaksPredicate implements IPredicate<ISpectrum>{

    public static final int MINIMUM_NUMBER_OF_PEAKS = 100;

    private final int miniNumOfPeaks;

    public MinimumNumberOfPeaksPredicate() {
        this(MINIMUM_NUMBER_OF_PEAKS);
    }

    public MinimumNumberOfPeaksPredicate(int miniNumOfPeaks) {
        this.miniNumOfPeaks = miniNumOfPeaks;
    }

    @Override
    public boolean apply(ISpectrum spectrum) {
        return spectrum.getPeaksCount() >= miniNumOfPeaks;
    }
}

package uk.ac.ebi.pride.spectracluster.util.predicate.spectrum_comparison;

import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.util.predicate.IComparisonPredicate;

import java.util.HashSet;
import java.util.Set;

/**
 * Created by jg on 20.05.15.
 */
public class ShareMajorPeaksPredicate implements IComparisonPredicate<ISpectrum> {
    public static final int DEFAULT_MAJOR_PEAKS = 5;
    private final int nMajorPeaks;

    public ShareMajorPeaksPredicate(int nMajorPeaks) {
        this.nMajorPeaks = nMajorPeaks;
    }

    public ShareMajorPeaksPredicate() {
        this(DEFAULT_MAJOR_PEAKS);
    }

    @Override
    public boolean apply(ISpectrum o1, ISpectrum o2) {
        int[] majorPeaks1 = o1.asMajorPeakMZs(nMajorPeaks);
        int[] majorPeaks2 = o2.asMajorPeakMZs(nMajorPeaks);

        Set<Integer> majorPeaks1Set = new HashSet<Integer>();
        for (int majorPeak : majorPeaks1)
            majorPeaks1Set.add(majorPeak);

        // check if any of the others exist
        for (int majorPeak : majorPeaks2) {
            if (majorPeaks1Set.contains(majorPeak)) {
                return true;
            }
        }

        return false;
    }
}

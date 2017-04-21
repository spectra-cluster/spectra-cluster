package uk.ac.ebi.pride.spectracluster.util.predicate.spectrum_comparison;

import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.util.predicate.IComparisonPredicate;

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

        // major peaks are always sorted according to m/z
        for (int index1 = 0, index2 = 0; index1 < majorPeaks1.length && index2 < majorPeaks2.length; ) {
            if (majorPeaks1[index1] == majorPeaks2[index2])
                return true;

            if (majorPeaks1[index1] < majorPeaks2[index2]) {
                index1++;
            }
            else  if (majorPeaks1[index1] > majorPeaks2[index2]) {
                index2++;
            }
            else {
                // both are identical, increment the lower one
                if (index1 < index2)
                    index1++;
                else
                    index2++;
            }
        }

        return false;
    }
}

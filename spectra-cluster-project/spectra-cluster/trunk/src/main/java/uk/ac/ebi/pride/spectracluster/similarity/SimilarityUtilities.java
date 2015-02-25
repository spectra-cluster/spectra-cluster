package uk.ac.ebi.pride.spectracluster.similarity;

import uk.ac.ebi.pride.spectracluster.spectrum.IPeak;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.spectrum.Peak;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by jg on 23.01.15.
 */
public final class SimilarityUtilities {
    private static final IPeak LAST_PEAK = new Peak(Float.MAX_VALUE, 0);

    protected SimilarityUtilities() {

    }

    public static PeakMatches getSharedPeaksAsMatches(ISpectrum spectrum1, ISpectrum spectrum2, float mzTolerance) {
        List<Integer>[] sharedPeakIndices = getSharedPeaks(spectrum1, spectrum2, mzTolerance);
        return new PeakMatches(spectrum1, spectrum2, sharedPeakIndices[0], sharedPeakIndices[1]);
    }

    /**
     * Finds the peaks shared between two spectra. This function returns the maximal number of
     * matches possible based on the set mzTolerance.
     * @param spectrum1 The first spectrum to match the peaks from.
     * @param spectrum2 The second spectrum to match the peaks from.
     * @param mzTolerance Peak tolerance for matching in m/z
     * @return Array of size 2. First list contains the peak indexes of spectrum 1, the second item the
     * corresponding indices of spectrum 2.
     */
    public static List<Integer>[] getSharedPeaks(ISpectrum spectrum1, ISpectrum spectrum2, float mzTolerance) {
        List<IPeak> peaks1 = spectrum1.getPeaks();
        List<IPeak> peaks2 = spectrum2.getPeaks();

        List<Integer> sharedPeaksIndexes1 = new ArrayList<Integer>();
        List<Integer> sharedPeaksIndexes2 = new ArrayList<Integer>();

        // upper and lower bound
        int indexSpec1 = 0, indexSpec2 = 0;

        while (indexSpec1 < peaks1.size() && indexSpec2 < peaks2.size()) {
            float mz1 = peaks1.get(indexSpec1).getMz();
            float mz2 = peaks2.get(indexSpec2).getMz();
            float difference = Math.abs(mz1 - mz2);

            if (difference > mzTolerance) {
                if (mz1 < mz2) {
                    indexSpec1++;
                    continue;
                }
                else {
                    indexSpec2++;
                    continue;
                }
            }
            // a potential match was found
            else {
                float differenceNextSpec1 = (indexSpec1 < peaks1.size() - 1) ?
                        Math.abs(peaks1.get(indexSpec1 + 1).getMz() - mz2) :
                        Float.MAX_VALUE;
                float differenceNextSpec2 = (indexSpec2 < peaks2.size() - 1) ?
                        Math.abs(peaks2.get(indexSpec2 + 1).getMz() - mz1) :
                        Float.MAX_VALUE;
                float differenceNextSpec1Spec2 = (indexSpec1 < peaks1.size() - 1 && indexSpec2 < peaks2.size() - 1) ?
                        Math.abs(peaks1.get(indexSpec1 + 1).getMz() - peaks2.get(indexSpec2 + 1).getMz()) :
                        Float.MAX_VALUE;

                // if the next two peaks are also a match, just match the current two
                if (differenceNextSpec1Spec2 <= mzTolerance) {
                    sharedPeaksIndexes1.add(indexSpec1);
                    sharedPeaksIndexes2.add(indexSpec2);

                    indexSpec1++;
                    indexSpec2++;
                    continue;
                }

                // using next peak in spec 1 is the best match
                if (differenceNextSpec1 < difference && differenceNextSpec1 < differenceNextSpec2) {
                    indexSpec1++;

                    sharedPeaksIndexes1.add(indexSpec1);
                    sharedPeaksIndexes2.add(indexSpec2);

                    indexSpec1++;
                    indexSpec2++;
                    continue;
                }

                // using next peak in spec 2 is the best match
                if (differenceNextSpec2 < difference) {
                    indexSpec2++;

                    sharedPeaksIndexes1.add(indexSpec1);
                    sharedPeaksIndexes2.add(indexSpec2);

                    indexSpec2++;
                    indexSpec1++;
                    continue;
                }

                // match the two current peaks
                sharedPeaksIndexes1.add(indexSpec1);
                sharedPeaksIndexes2.add(indexSpec2);

                indexSpec1++;
                indexSpec2++;
                continue;
            }
        }

        List<Integer>[] result = new List[2];
        result[0] = sharedPeaksIndexes1;
        result[1] = sharedPeaksIndexes2;

        return result;
    }

    /**
     * This function replicates the behaviour of the previous peak matching function used
     * in the FrankEtAlDotProduct class. It is focused on finding an optimal match. Thereby
     * peak matches within the tolerance may be lost.
     * @param spectrum1
     * @param spectrum2
     * @param mzTolerance
     * @return
     */
    @Deprecated
    public static List<Integer>[] getSharedPeaks2(ISpectrum spectrum1, ISpectrum spectrum2, float mzTolerance) {
        IPeak[] peaks1 = spectrum1.getPeaks().toArray(new IPeak[spectrum1.getPeaks().size() + 1]);
        IPeak[] peaks2 = spectrum2.getPeaks().toArray(new IPeak[spectrum2.getPeaks().size() + 1]);

        peaks1[peaks1.length - 1] = LAST_PEAK;
        peaks2[peaks2.length - 1] = LAST_PEAK;

        List<Integer> sharedPeaksSpec1 = new ArrayList<Integer>();
        List<Integer> sharedPeaksSpec2 = new ArrayList<Integer>();

        boolean spec1IsLast = false;
        int indexSpec1 = 0;
        int indexSpec2 = 0;

        while (indexSpec1 < peaks1.length - 1 && indexSpec2 < peaks2.length - 1) {
            IPeak peak1 = peaks1[indexSpec1];
            IPeak peak2 = peaks2[indexSpec2];

            double mz1 = peak1.getMz();
            double mz2 = peak2.getMz();
            double mzDifference = Math.abs(mz1 - mz2);

            if (Math.abs(mzDifference) <= mzTolerance) {
                // also calculate the difference for the next t and e peaks
                IPeak nextPeak1 = peaks1[indexSpec1 + 1];
                IPeak nextPeak2 = peaks2[indexSpec2 + 1];
                double mzDifferenceNextSpec1 = Math.abs(mz2 - nextPeak1.getMz());
                double mzDifferenceNextSpec2 = Math.abs(nextPeak2.getMz() - mz1);
                double mzDifferenceNextSpec12 = Math.abs(nextPeak2.getMz() - nextPeak1.getMz());

                // use the next peak in spectrum 1 if it's a better match
                if (mzDifferenceNextSpec1 < mzDifference &&
                    mzDifferenceNextSpec1 < mzDifferenceNextSpec2 &&
                    mzDifferenceNextSpec1 < mzDifferenceNextSpec12) {
                    indexSpec1++;
                }

                // use the next peak in spectrum 2 if it's a better match
                if (mzDifferenceNextSpec2 < mzDifference &&
                    mzDifferenceNextSpec2 < mzDifferenceNextSpec1 &&
                    mzDifferenceNextSpec2 < mzDifferenceNextSpec12) {
                    indexSpec2++;
                }

                sharedPeaksSpec1.add(indexSpec1);
                sharedPeaksSpec2.add(indexSpec2);

                // increment both counters since both peaks must not be compared again
                indexSpec1++;
                indexSpec2++;
                continue;
            }

            if (mzDifference == 0) {
                if (spec1IsLast) {
                    indexSpec2++;
                    spec1IsLast = false;
                } else {
                    indexSpec1++;
                    spec1IsLast = true;
                }
            } else {
                if (mz1 < mz2) {
                    indexSpec1++;
                } else {
                    indexSpec2++;
                }
            }

        }

        List<Integer>[] result = new List[2];
        result[0] = sharedPeaksSpec1;
        result[1] = sharedPeaksSpec2;

        return result;
    }
}

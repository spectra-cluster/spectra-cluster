package uk.ac.ebi.pride.spectracluster.similarity;

import uk.ac.ebi.pride.spectracluster.spectrum.IPeak;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;

import java.util.List;

/**
 * This class holds information about the matched peaks between
 * two spectra. It is intended to be used if multiple
 * ISimilarityCheckerS are used on the same spectrum pairs.
 * Thereby, the peak matching process only has to be done ones.
 *
 * Created by jg on 25.02.15.
 */
public class PeakMatches {
    private final ISpectrum spectrum1;
    private final ISpectrum spectrum2;
    private final List<Integer> sharedPeakIndecesSpec1;
    private final List<Integer> sharedPeakIndecesSpec2;

    private IPeak[] sharedPeaksSpec1 = null;
    private IPeak[] sharedPeaksSpec2 = null;

    public PeakMatches(ISpectrum spectrum1, ISpectrum spectrum2, List<Integer> sharedPeakIndecesSpec1, List<Integer> sharedPeakIndecesSpec2) {
        this.spectrum1 = spectrum1;
        this.spectrum2 = spectrum2;
        this.sharedPeakIndecesSpec1 = sharedPeakIndecesSpec1;
        this.sharedPeakIndecesSpec2 = sharedPeakIndecesSpec2;

        // make sure the number of indices have the same length
        if (sharedPeakIndecesSpec1.size() != sharedPeakIndecesSpec2.size())
            throw new IllegalStateException("SharedPeakIndices must be of same size for both spectra");
    }

    /**
     * Returns the shared peaks from spectrum 1
     * @return
     */
    public IPeak[] getSharedPeaksSpec1() {
        if (sharedPeaksSpec1 == null) {
            sharedPeaksSpec1 = new IPeak[sharedPeakIndecesSpec1.size()];

            for (int i = 0; i < sharedPeakIndecesSpec1.size(); i++) {
                sharedPeaksSpec1[i] = spectrum1.getPeaks().get(sharedPeakIndecesSpec1.get(i));
            }
        }

        return sharedPeaksSpec1;
    }

    /**
     * Returns the shared peaks from spectrum 2
     * @return
     */
    public IPeak[] getSharedPeaksSpec2() {
        if (sharedPeaksSpec2 == null) {
            sharedPeaksSpec2 = new IPeak[sharedPeakIndecesSpec2.size()];

            for (int i = 0; i < sharedPeakIndecesSpec2.size(); i++) {
                sharedPeaksSpec2[i] = spectrum2.getPeaks().get(sharedPeakIndecesSpec2.get(i));
            }
        }

        return sharedPeaksSpec2;
    }

    /**
     * Returns the number of shared peaks
     * @return
     */
    public int getNumberOfSharedPeaks() {
        return sharedPeakIndecesSpec1.size();
    }

    /**
     * Returns the defined peak pair (0-based index). Use
     * 'getNumberOfSharedPeaks' to determine the number of
     * shared peaks.
     * @param nIndex 0-based index of the peak pair.
     * @return IPeak[] - index 0: peak from spectrum 1, index1: peak from spectrum 2
     */
    public IPeak[] getPeakPair(int nIndex) {
        if (nIndex < 0)
            throw new IndexOutOfBoundsException("PeakPair index must be greater than 0");
        if (nIndex >= sharedPeakIndecesSpec1.size())
            throw new IndexOutOfBoundsException("Request PeakPair with index '" + nIndex + "' from " + sharedPeakIndecesSpec1.size() + " matches");

        IPeak[] peakPair = new IPeak[2];
        peakPair[0] = spectrum1.getPeaks().get(sharedPeakIndecesSpec1.get(nIndex));
        peakPair[1] = spectrum2.getPeaks().get(sharedPeakIndecesSpec2.get(nIndex));

        return peakPair;
    }

    public ISpectrum getSpectrum1() {
        return spectrum1;
    }

    public ISpectrum getSpectrum2() {
        return spectrum2;
    }
}

package uk.ac.ebi.pride.spectracluster.similarity;

import uk.ac.ebi.pride.spectracluster.spectrum.IPeak;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.util.Pair;

import java.util.List;

/**
 * This interface holds information about the matched peaks between
 * two spectra. It is intended to be used if multiple
 * ISimilarityCheckers are used on the same spectrum pairs.
 * Thereby, the peak matching process only has to be done ones.
 *
 * @author Johannes Griss
 * @author Rui Wang
 * @version $Id$
 */
public interface IPeakMatches {

    /**
     * Get shared peaks from spectrum one
     *
     * @return an array of peaks
     */
    List<IPeak> getSharedPeaksFromSpectrumOne();

    /**
     * Get shared peaks from spectrum two
     *
     * @return an array of peaks
     */
    List<IPeak> getSharedPeaksFromSpectrumTwo();

    /**
     * Get the number of shared peaks
     *
     * @return the number of shared peaks
     */
    int getNumberOfSharedPeaks();

    /**
     * Get a pair of peaks at a given index
     *
     * @param nIndex index of the shared peak
     * @return an array of size two containing the peaks
     */
    Pair<IPeak, IPeak> getPeakPair(int nIndex);

    /**
     * Get spectrum one
     *
     * @return spectrum
     */
    ISpectrum getSpectrumOne();

    /**
     * Get spectrum two
     *
     * @return spectrum
     */
    ISpectrum getSpectrumTwo();
}

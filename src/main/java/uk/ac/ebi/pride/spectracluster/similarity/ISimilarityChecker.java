package uk.ac.ebi.pride.spectracluster.similarity;

import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.util.IAlgorithm;

/**
 * Assesses the similarity between two
 * spectra and returns the result of this
 * check as a double. Higher values mean
 * a higher similarity.
 *
 * @author jg
 * @author Rui Wang
 */
public interface ISimilarityChecker extends IAlgorithm {

    /**
     * Assesses the similarity between the two passed
     * spectra and returns the result of this assessment
     * as a double. A higher number reflects a higher
     * grade of similarity.
     *
     * @param spectrum1 The first spectrum to compare. The list of Peaks MUST be sorted according to intensity.
     * @param spectrum2 The second spectrum to compare. The list of Peaks MUST be sorted according to intensity.
     * @return A score indicating the similarity between the two passed spectra.
     */
    double assessSimilarity(ISpectrum spectrum1, ISpectrum spectrum2);

    double assessSimilarity(IPeakMatches peakMatches);

    /**
     * Indicates whether peak filtering is enabled for the
     * current algorithm
     * @return
     */
    boolean isPeakFiltering();

    /**
     * Allows the disabling of peak filtering for the current
     * algorithm.
     * @param peakFiltering
     */
    void setPeakFiltering(boolean peakFiltering);

    /**
     * Set the used fragment ion tolerance used for peak matching between the two spectra.
     * @param fragmentIonTolerance Fragment tolerance in m/z
     */
    void setFragmentIonTolerance(float fragmentIonTolerance);

    /**
     * Returns the currently used fragment ion tolerance.
     * @return The fragment ion tolerance in m/z
     */
    float getFragmentIonTolerance();
}

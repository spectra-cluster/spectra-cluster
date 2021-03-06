package uk.ac.ebi.pride.tools.pride_spectra_clustering.normalizer;

import uk.ac.ebi.pride.tools.pride_spectra_clustering.util.Peak;

import java.util.List;

/**
 * Normalizes a spectrum's intensities.
 *
 * @author jg
 */
public interface IntensityNormalizer {
    /**
     * Normalizes the given spectrum's intensities.
     *
     * @param spectrum The spectrum as a Map with the m/z values as key and their intensities as values.
     * @return The normalized version of the spectrum as a map.
     */
    List<Peak> normalizeSpectrum(List<Peak> spectrum);
}

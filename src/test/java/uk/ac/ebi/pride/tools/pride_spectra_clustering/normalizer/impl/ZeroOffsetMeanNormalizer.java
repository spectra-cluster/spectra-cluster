package uk.ac.ebi.pride.tools.pride_spectra_clustering.normalizer.impl;

import uk.ac.ebi.pride.tools.pride_spectra_clustering.normalizer.IntensityNormalizer;
import uk.ac.ebi.pride.tools.pride_spectra_clustering.util.Peak;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Normalizes a spectrum by dividing the intensities
 * through the mean of all intensities in the spectrum.
 * TODO: add reference
 *
 * @author jg
 */
public class ZeroOffsetMeanNormalizer implements IntensityNormalizer {

    public List<Peak> normalizeSpectrum(List<Peak> spectrum) {
        // calculate the total intensity
        Double totalIntensity = 0.0;

        for (Peak p : spectrum)
            totalIntensity += p.getIntensity();

        double averageIntensity = totalIntensity / spectrum.size();

        // normalize the spectrum
        List<Peak> normalizedSpectrum = spectrum.stream()
                .map(p -> new Peak(p.getMz(), p.getIntensity() / averageIntensity))
                .collect(Collectors.toCollection(() -> new ArrayList<>(spectrum.size())));

        return normalizedSpectrum;
    }

}

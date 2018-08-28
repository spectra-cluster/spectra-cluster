package uk.ac.ebi.pride.tools.pride_spectra_clustering.normalizer.impl;

import uk.ac.ebi.pride.tools.pride_spectra_clustering.normalizer.IntensityNormalizer;
import uk.ac.ebi.pride.tools.pride_spectra_clustering.util.Peak;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Normalizes a spectrum by dividing
 * the intensities through the square
 * root of the sum of all squared
 * intensities.
 *
 * @author jg
 */
public class SquareRootSumNormalizer implements IntensityNormalizer {

    public List<Peak> normalizeSpectrum(List<Peak> spectrum) {
        // get the intensities
        List<Double> intensities = spectrum.stream()
                .map(Peak::getIntensity)
                .collect(Collectors.toCollection(() -> new ArrayList<>(spectrum.size())));

        // calculate the sum of the squared intensities
        double sumSquaredIntensities = intensities.stream()
                .mapToDouble(intensity -> Math.pow(intensity, 2))
                .sum();

        // calculate the ratio by taking the square root of the summed intensities
        double ratio = Math.sqrt(sumSquaredIntensities);

        // create the new spectrum
        List<Peak> normalizedSpectrum = spectrum.stream()
                .map(p -> new Peak(p.getMz(), p.getIntensity() * ratio))
                .collect(Collectors.toCollection(() -> new ArrayList<>(spectrum.size())));

        return normalizedSpectrum;
    }

}

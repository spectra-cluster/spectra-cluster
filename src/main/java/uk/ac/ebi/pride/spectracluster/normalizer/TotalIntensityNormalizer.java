package uk.ac.ebi.pride.spectracluster.normalizer;

import uk.ac.ebi.pride.spectracluster.spectrum.IPeak;
import uk.ac.ebi.pride.spectracluster.spectrum.Peak;
import uk.ac.ebi.pride.spectracluster.util.CompareTo;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Normalizes a spectrum's intensities so that
 * the spectrum's total intensity matches a given
 * number.
 * This method was used by Frank et al. (2008) JPR: Clustering millions of spectra
 * To calculate the dot-product they then used
 * 1+ln(intensity) as the peak's intensity
 *
 * @author jg
 * @author Rui Wang
 */
public class TotalIntensityNormalizer implements IIntensityNormalizer {
    private static final String VERSION = "1.0";

    private static final double DEFAULT_TOTAL_INTENSITY = 1000;

    private final double totalIntensity;

    public static final IIntensityNormalizer DEFAULT = new TotalIntensityNormalizer();

    /**
     * Use Defaults which builds with reflection
     * Set the class with Defaults.setNormalizerClass
     * use DEFAULT
     */
    public TotalIntensityNormalizer() {
        this(DEFAULT_TOTAL_INTENSITY);
    }

    /**
     * Use Defaults which builds with reflection
     * Set the class with Defaults.setNormalizerClass
     */
    public TotalIntensityNormalizer(double total) {
        totalIntensity = total;
    }

    /**
     * return a name which should not change
     *
     * @return !null name
     */
    @Override
    public String getName() {
        return getClass().getSimpleName();
    }

    /**
     * return a version number - this may be updated over time
     *
     * @return !null version
     */
    @Override
    public String getCurrentVersion() {
        return VERSION;
    }


    public double getTotalIntensity() {
        return totalIntensity;
    }

    /**
     * return the value normalized to - especial;ly useful for total intensity normalization where
     * we may not weed to normalize
     *
     * @return as above
     */
    @Override
    public double getNormalizedValue() {
        return getTotalIntensity();
    }

    @Override
    public List<IPeak> normalizePeaks(List<IPeak> peaks) {
        // get the max intensity
        Double specTotalIntensity = 0.0;

        // create the new spectrum
        List<IPeak> normalizedSpectrum = new ArrayList<>(peaks.size());


        for (IPeak p : peaks) {
            specTotalIntensity += p.getIntensity();
        }

        // if there's no suitable max intensity, return the unchanged spectrum
        if (specTotalIntensity <= 0)
            return normalizedSpectrum;

        // calculate the ratio
        double ratio = getTotalIntensity() / specTotalIntensity;


        normalizedSpectrum = peaks.stream()
                .map(p -> new Peak(p.getMz(), (float) (p.getIntensity() * ratio)))
                .collect(Collectors.toCollection(() -> new ArrayList<>(peaks.size())));

        return normalizedSpectrum;
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        final TotalIntensityNormalizer that = (TotalIntensityNormalizer) o;

        return CompareTo.compare(that.totalIntensity, totalIntensity) == 0;

    }

    @Override
    public int hashCode() {
        final long temp = Double.doubleToLongBits(totalIntensity);
        return (int) (temp ^ (temp >>> 32));
    }
}

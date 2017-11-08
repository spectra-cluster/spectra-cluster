package uk.ac.ebi.pride.spectracluster.quality;

import uk.ac.ebi.pride.spectracluster.spectrum.IPeak;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;

import java.util.ArrayList;
import java.util.Collections;

/**
 * Assesses a spectrum's signal-to-noise
 * ratio as a measure of the spectrum's
 * quality. The method used is taken from
 * Lam et al. (2008), Nat Methods 5(10):873-875
 *
 * @author jg
 */
public class OriginalSignalToNoiseChecker implements IQualityScorer {
    public static final String VERSION = "1.0";

    public static final int NUMBER_SPECTRA_TO_CHECK = 6;


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

    /**
     * Calculates a spectrum's signal-to-noise ratio
     * by taking the 2nd-6th highest peak's intensity
     * and dividing it by the median intensity of all
     * peaks.
     * This method should only be used on normalized
     * spectra.
     */
    @Override
    public double calculateQualityScore(ISpectrum spectrum) {
        // get the intensities
        ArrayList<Float> intensities = new ArrayList<>(spectrum.getPeaks().size());
        for (IPeak p : spectrum.getPeaks())
            intensities.add(p.getIntensity());

        Collections.sort(intensities);

        // get the total intensity of the 2nd-6th highest peak
        double highestPeakIntensity = 0.0;

        if (intensities.size() < NUMBER_SPECTRA_TO_CHECK)
            return 0.0;

        for (int n = intensities.size() - NUMBER_SPECTRA_TO_CHECK; n < intensities.size() - 1; n++)
            highestPeakIntensity += intensities.get(n);

        highestPeakIntensity = highestPeakIntensity / (NUMBER_SPECTRA_TO_CHECK - 1);

        // get the median
        double nPeaks = intensities.size();
        double median = 0.0;

        // check if there's an even number of peaks
        if (nPeaks % 2 == 0) {
            int index1 = (int) nPeaks / 2 - 1;
            Float intensity1 = intensities.get(index1);
            int index2 = (int) nPeaks / 2;
            Float intensity2 = intensities.get(index2);
            median = (intensity1 + intensity2) / 2;
        } else {
            int index = (int) nPeaks / 2;
            median = intensities.get(index);
        }

        return highestPeakIntensity / median;
    }
}

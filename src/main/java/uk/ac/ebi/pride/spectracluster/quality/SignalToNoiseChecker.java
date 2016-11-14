package uk.ac.ebi.pride.spectracluster.quality;

import uk.ac.ebi.pride.spectracluster.spectrum.IPeak;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.util.comparator.PeakIntensityComparator;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * uk.ac.ebi.pride.spectracluster.quality.SignalToNoiseChecker
 *
 * @author Steve Lewis
 */
public class SignalToNoiseChecker implements IQualityScorer {

    private static final int NUMBER_HIGH_PEAKS = 6;
    private static final String VERSION = "1.0";


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
        ISpectrum highestNPeaks = spectrum.getHighestNPeaks(NUMBER_HIGH_PEAKS);

        if (highestNPeaks.getPeaksCount() < NUMBER_HIGH_PEAKS)
            return 0.0;
        double totalIntensity = highestNPeaks.getTotalIntensity();
        double highestPeak = 0;
        for (IPeak peak : highestNPeaks.getPeaks()) {
            highestPeak = Math.max(peak.getIntensity(), highestPeak);
        }

        double meanHigh = (totalIntensity - highestPeak) / (NUMBER_HIGH_PEAKS - 1);

        List<IPeak> peaks = new ArrayList<IPeak>(spectrum.getPeaks());
        Collections.sort(peaks, PeakIntensityComparator.INSTANCE);

        double median;

        int peakSize = peaks.size();
        if (peakSize % 2 == 1) {
            int index = peakSize / 2;
            IPeak iPeak = peaks.get(index);
            median = iPeak.getIntensity(); // odd case

        } else {
            int index2 = (peakSize / 2);
            int index1 = index2 - 1;
            IPeak iPeak1 = peaks.get(index1);
            IPeak iPeak2 = peaks.get(index2);
            double intensity1 = iPeak1.getIntensity();
            double intensity2 = iPeak2.getIntensity();
            median = (intensity1 + intensity2) / 2; // even case
        }

        return meanHigh / median;
    }
}

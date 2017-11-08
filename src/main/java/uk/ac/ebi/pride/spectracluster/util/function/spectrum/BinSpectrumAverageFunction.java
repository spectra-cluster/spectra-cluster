package uk.ac.ebi.pride.spectracluster.util.function.spectrum;

import uk.ac.ebi.pride.spectracluster.spectrum.IPeak;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.spectrum.Peak;
import uk.ac.ebi.pride.spectracluster.spectrum.Spectrum;
import uk.ac.ebi.pride.spectracluster.util.function.IFunction;

import java.util.ArrayList;
import java.util.List;

/**
 * Bins the spectrum's peaks and merges all peaks
 * into a single one based on the set tolerance.
 *
 * When merging peaks the new m/z is weighted based on the
 * peaks' intensities. The intensities are summed and the
 * maximum count of the two peak is used.
 *
 * This function requires the peaks to be sorted by m/z.
 *
 * Created by jg on 06.11.17.
 */
public class BinSpectrumAverageFunction implements IFunction<ISpectrum, ISpectrum> {
    private final float fragmentTolerance;

    /**
     * Creates a new BinSpectrumFunction.
     * @param fragmentTolerance The tolerance to use to merge peaks.
     */
    public BinSpectrumAverageFunction(float fragmentTolerance) {
        this.fragmentTolerance = fragmentTolerance;
    }

    @Override
    public ISpectrum apply(ISpectrum o) {
        List<IPeak> binnedPeaks = new ArrayList<IPeak>();
        IPeak currentPeak = null;

        for (IPeak peak : o.getPeaks()) {
            if (currentPeak == null) {
                currentPeak = peak;
                continue;
            }

            // make sure the sort order is correct - this should never happen
            if (peak.getMz() < currentPeak.getMz()) {
                throw new IllegalStateException("BinSpectrumFunction can only be applied to spectra with sorted peaks.");
            }

            // merge peaks within the tolerance
            if (peak.getMz() < currentPeak.getMz() + fragmentTolerance) {
                float totalIntensity = peak.getIntensity() + currentPeak.getIntensity();
                float weightedMz = (peak.getMz() * peak.getIntensity() / totalIntensity) +
                        currentPeak.getMz() * currentPeak.getIntensity() / totalIntensity;
                // use the max count that was observed until now
                int maxCount = Math.max(currentPeak.getCount(), peak.getCount());

                currentPeak = new Peak(weightedMz, currentPeak.getIntensity() + peak.getIntensity(), maxCount);
            } else {
                // create a new peak
                binnedPeaks.add(currentPeak);
                currentPeak = peak;
            }
        }

        if (currentPeak != null) {
            binnedPeaks.add(currentPeak);
        }

        // create the new spectrum
        return new Spectrum(o, binnedPeaks, true);
    }
}

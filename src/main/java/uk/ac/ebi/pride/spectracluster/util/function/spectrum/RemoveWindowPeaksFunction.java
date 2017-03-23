package uk.ac.ebi.pride.spectracluster.util.function.spectrum;

import uk.ac.ebi.pride.spectracluster.spectrum.IPeak;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.spectrum.Masses;
import uk.ac.ebi.pride.spectracluster.spectrum.Spectrum;
import uk.ac.ebi.pride.spectracluster.util.function.IFunction;

import java.util.ArrayList;
import java.util.List;

/**
 * This filter removes all peaks that are
 * outside of a defined window. The default
 * version removes all peaks below 200 m/z.
 *
 * Created by jg on 13.05.15.
 */
public class RemoveWindowPeaksFunction implements IFunction<ISpectrum, ISpectrum> {
    private final float minMz;
    private final float maxMz;

    /**
     * Removes all peaks that are outside the defined
     * limits.
     * @param minMz The minimum m/z peaks must have.
     * @param maxMz The maximum m/z peaks may have.
     */
    public RemoveWindowPeaksFunction(float minMz, float maxMz) {
        this.minMz = minMz;
        this.maxMz = maxMz;
    }

    /**
     * Default constructor removes all peaks
     * below 200 m/z.
     */
    public RemoveWindowPeaksFunction() {
        this.minMz = 200F;
        this.maxMz = Float.MAX_VALUE;
    }

    @Override
    public ISpectrum apply(ISpectrum o) {
        List<IPeak> filteredPeaks = new ArrayList<IPeak>();

        for (IPeak peak : o.getPeaks()) {
            if (peak.getMz() < minMz || peak.getMz() > maxMz)
                continue;

            filteredPeaks.add(peak);
        }

        ISpectrum filteredSpectrum = new Spectrum(o, filteredPeaks, true); // does not require resorting

        return filteredSpectrum;
    }
}

package uk.ac.ebi.pride.spectracluster.util.function.spectrum;

import uk.ac.ebi.pride.spectracluster.spectrum.IPeak;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.spectrum.Spectrum;
import uk.ac.ebi.pride.spectracluster.util.function.IFunction;
import uk.ac.ebi.pride.spectracluster.util.function.peak.HighestNPeakFunction;

import java.util.List;

/**
 * Created by jg on 13.05.15.
 */
public class HighestNSpectrumPeaksFunction implements IFunction<ISpectrum, ISpectrum> {
    public final IFunction<List<IPeak>, List<IPeak>> peakFilter;

    public HighestNSpectrumPeaksFunction(int maxPeaks) {
        peakFilter = new HighestNPeakFunction(maxPeaks);
    }

    @Override
    public ISpectrum apply(ISpectrum o) {
        return new Spectrum(o, peakFilter.apply(o.getPeaks()));
    }
}

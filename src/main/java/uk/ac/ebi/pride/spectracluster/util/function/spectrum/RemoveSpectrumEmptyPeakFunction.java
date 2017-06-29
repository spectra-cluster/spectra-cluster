package uk.ac.ebi.pride.spectracluster.util.function.spectrum;

import uk.ac.ebi.pride.spectracluster.spectrum.IPeak;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.spectrum.Spectrum;
import uk.ac.ebi.pride.spectracluster.util.function.IFunction;

import java.util.List;
import java.util.stream.Collectors;

/**
 * Remove peaks where intensities are zero
 *
 * @author Rui Wang
 * @version $Id$
 */
public class RemoveSpectrumEmptyPeakFunction implements IFunction<ISpectrum, ISpectrum> {

    @Override
    public ISpectrum apply(ISpectrum spectrum) {
        List<IPeak> peaks = spectrum.getPeaks();

        // vet each peak
        List<IPeak> vettedPeaks = peaks.stream()
                .filter(peak -> peak.getIntensity() > 0 && peak.getMz() > 0)
                .collect(Collectors.toList());

        return new Spectrum(spectrum, vettedPeaks, true);
    }
}

package uk.ac.ebi.pride.spectracluster.util.function.spectrum;

import uk.ac.ebi.pride.spectracluster.spectrum.IPeak;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.spectrum.Spectrum;
import uk.ac.ebi.pride.spectracluster.util.function.IFunction;

import java.util.ArrayList;
import java.util.List;

/**
 * Remove peaks where intensities are zero
 *
 * @author Rui Wang
 * @version $Id$
 */
public class RemoveEmptyPeakFunction implements IFunction<ISpectrum, ISpectrum> {

    @Override
    public ISpectrum apply(ISpectrum spectrum) {
        List<IPeak> peaks = spectrum.getPeaks();

        // vet each peak
        List<IPeak> vettedPeaks = new ArrayList<IPeak>();
        for (IPeak peak : peaks) {
            if (peak.getIntensity() > 0 && peak.getMz() > 0)
                vettedPeaks.add(peak);
        }

        return new Spectrum(spectrum, vettedPeaks);
    }
}

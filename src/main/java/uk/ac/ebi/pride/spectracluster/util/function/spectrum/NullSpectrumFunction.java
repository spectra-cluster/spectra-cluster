package uk.ac.ebi.pride.spectracluster.util.function.spectrum;

import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.spectrum.Spectrum;
import uk.ac.ebi.pride.spectracluster.util.function.IFunction;

/**
 * Created by jg on 13.05.15.
 */
public class NullSpectrumFunction implements IFunction<ISpectrum, ISpectrum> {
    @Override
    public ISpectrum apply(ISpectrum o) {
        return new Spectrum(o);
    }
}

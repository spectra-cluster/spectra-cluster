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
 * above the precursor mass + a tolerance
 *
 * Created by jg on 13.05.15.
 */
public class RemoveImpossiblyHighPeaksFunction implements IFunction<ISpectrum, ISpectrum> {
    public final static float DEFAULT_TOLERANCE = 3.0F;
    public final float tolerance;

    public RemoveImpossiblyHighPeaksFunction(float tolerance) {
        this.tolerance = tolerance;
    }

    public RemoveImpossiblyHighPeaksFunction() {
        this(DEFAULT_TOLERANCE);
    }

    @Override
    public ISpectrum apply(ISpectrum o) {
        final float monoisotopicMass = Masses.getMonoisotopicMass(o.getPrecursorMz(), o.getPrecursorCharge());
        final float maxMass = monoisotopicMass + Masses.PROTON + tolerance;

        List<IPeak> filteredPeaks = new ArrayList<IPeak>();

        for (IPeak peak : o.getPeaks()) {
            if (peak.getMz() > maxMass)
                continue;

            filteredPeaks.add(peak);
        }

        ISpectrum filteredSpectrum = new Spectrum(o, filteredPeaks);

        return filteredSpectrum;
    }
}

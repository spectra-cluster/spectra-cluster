package uk.ac.ebi.pride.spectracluster.util.function.spectrum;

import uk.ac.ebi.pride.spectracluster.spectrum.IPeak;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.spectrum.Masses;
import uk.ac.ebi.pride.spectracluster.spectrum.Spectrum;
import uk.ac.ebi.pride.spectracluster.util.function.IFunction;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by jg on 13.05.15.
 */
public class RemovePrecursorPeaksFunction implements IFunction<ISpectrum, ISpectrum> {
    private final float fragmentIonTolerance;

    public RemovePrecursorPeaksFunction(float fragmentIonTolerance) {
        this.fragmentIonTolerance = fragmentIonTolerance;
    }

    @Override
    public ISpectrum apply(ISpectrum o) {

        // calculate m/z of neutral losses
        final float floatCharge     = (float) o.getPrecursorCharge();
        final float waterLoss       = o.getPrecursorMz() - (Masses.WATER_MONO / floatCharge);
        final float doubleWaterLoss = o.getPrecursorMz() - (2.0F * Masses.WATER_MONO / floatCharge);
        final float ammoniumLoss    = o.getPrecursorMz() - (Masses.AMMONIA_MONO / floatCharge);

        // calculate range based on fragmentIonTolerance
        final float minWaterLoss        = waterLoss - fragmentIonTolerance;
        final float maxWaterLoss        = waterLoss + fragmentIonTolerance;
        final float minDoubleWaterLoss  = doubleWaterLoss - fragmentIonTolerance;
        final float maxDoubleWaterLoss  = doubleWaterLoss + fragmentIonTolerance;
        final float minAmmoniumLoss     = ammoniumLoss - fragmentIonTolerance;
        final float maxAmmoniumLoss     = ammoniumLoss + fragmentIonTolerance;

        // also filter the default precursor
        final float minPrecursor = o.getPrecursorMz() - fragmentIonTolerance;
        final float maxPrecursor = o.getPrecursorMz() + fragmentIonTolerance;

        // filter the peak list
        List<IPeak> filteredPeakList = new ArrayList<IPeak>();

        for (IPeak peak : o.getPeaks()) {
            final float peakMz = peak.getMz();

            // ignore any peak that could be a neutral loss
            if (isWithinRange(minWaterLoss, maxWaterLoss, peakMz))
                continue;
            if (isWithinRange(minDoubleWaterLoss, maxDoubleWaterLoss, peakMz))
                continue;
            if (isWithinRange(minAmmoniumLoss, maxAmmoniumLoss, peakMz))
                continue;
            if (isWithinRange(minPrecursor, maxPrecursor, peakMz))
                continue;

            filteredPeakList.add(peak);
        }

        ISpectrum filteredSpectrum = new Spectrum(o, filteredPeakList, true);

        return filteredSpectrum;
    }

    private boolean isWithinRange(float min, float max, float value) {
        return (value >= min && value <= max);
    }
}

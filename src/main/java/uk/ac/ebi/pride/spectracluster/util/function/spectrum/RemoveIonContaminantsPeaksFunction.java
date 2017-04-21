package uk.ac.ebi.pride.spectracluster.util.function.spectrum;

import uk.ac.ebi.pride.spectracluster.spectrum.IPeak;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.spectrum.Spectrum;
import uk.ac.ebi.pride.spectracluster.util.function.IFunction;

import java.util.ArrayList;
import java.util.List;

/**
 * Removes common ion contaminants from the spectrum
 * such as immonium ions.
 */
public class RemoveIonContaminantsPeaksFunction implements IFunction<ISpectrum, ISpectrum> {
    private final float fragmentIonTolerance;

    /**
     * For this filter function, a minimum fragment
     * ion tolerance of 0.1 m/z is recommended - mainly
     * because for many ions the monoisotopic mass was
     * not found
     * @param fragmentIonTolerance
     */
    public RemoveIonContaminantsPeaksFunction(float fragmentIonTolerance) {
        this.fragmentIonTolerance = fragmentIonTolerance;
    }

    public final float[] CONTAMINANT_ION_MZ = {
            /**
             * Immonium ions
             * This list of the most commonly observed
             * immonium ions is taken from this
             * page (accessed 23.02.2017):
             * http://www.ionsource.com/Card/immon/immon.htm
             * http://www.ionsource.com/Card/immon/more.htm
             */
            60.04494F, // Serine
            70.06568F, // Proline
            72.08133F, // Valine
            74.06059F, // Threonine
            86.09698F, // Leucine / Isoleucine
            87.05584F, // Asparagine
            88.03986F, // Aspartic acid
            101.0715F, // Glutamine
            102.0555F, // Glutamic acid
            104.0534F, // Methionine
            110.0718F, // Histidine
            120.0483F, // Ox Methionine
            120.0813F, // Phenylalanine
            133.0436F, // Carbamidomethylated C
            134.0276F, // Carboxymethylated C
            136.0762F, // Tyrosine
            147.0772F, // C + Acrylamide
            159.0922F, // Tryptophane

    };

    @Override
    public ISpectrum apply(ISpectrum o) {
        // filter the peak list
        List<IPeak> filteredPeakList = new ArrayList<IPeak>();

        for (IPeak peak : o.getPeaks()) {
            final float peakMz = peak.getMz();

            // ignore any peak that could be a neutral loss
            boolean isContaminantPeak = false;
            for (float contaminantMz : CONTAMINANT_ION_MZ) {
                if (isWithinRange(contaminantMz - fragmentIonTolerance, contaminantMz + fragmentIonTolerance, peakMz)) {
                    isContaminantPeak = true;
                    break;
                }
            }
            if (!isContaminantPeak) {
                filteredPeakList.add(peak);
            }
        }

        return new Spectrum(o, filteredPeakList, true);
    }

    private boolean isWithinRange(float min, float max, float value) {
        return (value >= min && value <= max);
    }
}

package uk.ac.ebi.pride.spectracluster.util.comparator;

import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;

import java.util.Comparator;

/**
 * Compare spectra by ID
 *
 * @author Rui Wang
 * @version $Id$
 */
public class SpectrumIDComparator implements Comparator<ISpectrum> {

    public static SpectrumIDComparator INSTANCE = new SpectrumIDComparator();

    private SpectrumIDComparator() {
    }

    @Override
    public int compare(ISpectrum o1, ISpectrum o2) {
        return o1.getId().compareTo(o2.getId());
    }
}

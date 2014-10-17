package uk.ac.ebi.pride.spectracluster.util.comparator;

import uk.ac.ebi.pride.spectracluster.spectrum.IPeak;

import java.util.Comparator;

/**
 * Comparator to compare peaks by intensity first then mz rather than the
 * standard mz then intensity
 *
 * @author Steve Lewis
 * @author Rui Wang
 * @version $Id$
 */
public final class PeakIntensityComparator implements Comparator<IPeak> {

    public static PeakIntensityComparator INSTANCE = new PeakIntensityComparator();

    private PeakIntensityComparator() {
    }

    /**
     * this version of the comparator puts the highest intensity peaks at the top
     *
     * @param o1
     * @param o2
     * @return
     */
    public int compare(IPeak o1, IPeak o2) {
        if (o1 == null) {
            return o2 == null ? 0 : -1;
        }

        if (o2 == null) {
            return 1;
        }

        if (o1.getIntensity() != o2.getIntensity()) {
            return o2.getIntensity() < o1.getIntensity() ? -1 : 1;
        }

        return 0;
    }
}

package uk.ac.ebi.pride.spectracluster.util.comparator;

import uk.ac.ebi.pride.spectracluster.spectrum.IPeak;

import java.util.Comparator;

/**
 * @author Rui Wang
 * @version $Id$
 */
public class PeakMzComparator implements Comparator<IPeak> {

    public int compare(IPeak o1, IPeak o2) {
        // check whether one of the peaks == null
        if (o1 == null && o2 == null) {
            return 0;
        }

        if (o1 == null) {
            return -1;
        }

        if (o2 == null) {
            return 1;
        }

        // sort according to m/z
        return Float.compare(o1.getMz(), o2.getMz());
    }
}
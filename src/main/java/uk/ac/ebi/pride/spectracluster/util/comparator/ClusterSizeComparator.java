package uk.ac.ebi.pride.spectracluster.util.comparator;

import uk.ac.ebi.pride.spectracluster.cluster.ICluster;

import java.util.Comparator;

/**
 * Compare clusters by the number of spectra within
 *
 * @author Rui Wang
 * @version $Id$
 */
public class ClusterSizeComparator implements Comparator<ICluster> {

    @Override
    public int compare(ICluster o1, ICluster o2) {
        int diff = o1.getClusteredSpectraCount() - o2.getClusteredSpectraCount();
        if (diff != 0)
            return diff > 0 ? -1 : 1;
        return o1.compareTo(o2);
    }
}

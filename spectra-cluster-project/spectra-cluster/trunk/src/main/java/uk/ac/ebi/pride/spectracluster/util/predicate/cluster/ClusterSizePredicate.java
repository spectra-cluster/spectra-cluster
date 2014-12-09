package uk.ac.ebi.pride.spectracluster.util.predicate.cluster;

import uk.ac.ebi.pride.spectracluster.cluster.ICluster;
import uk.ac.ebi.pride.spectracluster.util.predicate.IPredicate;

/**
 * Assess the size of the cluster. True if the number of spectra is more the specified value
 *
 * @author Rui Wang
 * @version $Id$
 */
public class ClusterSizePredicate implements IPredicate<ICluster> {

    public static final int DEFAULT_MINIMUM_CLUSTER_SIZE = 10;

    private final int miniSize;

    public ClusterSizePredicate() {
        this(DEFAULT_MINIMUM_CLUSTER_SIZE);
    }

    public ClusterSizePredicate(int miniSize) {
        this.miniSize = miniSize;
    }

    @Override
    public boolean apply(ICluster cluster) {
        return cluster.getClusteredSpectraCount() >= miniSize;
    }
}

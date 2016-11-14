package uk.ac.ebi.pride.spectracluster.cluster;

/**
 * IClusterStabilityAssessor defines an interface for assessing the stability of a cluster
 *
 * @author Rui Wang
 * @version $Id$
 */
public interface
        IClusterStabilityAssessor {

    /**
     * if true the cluster is stable and will not allow removal
     *
     * @return Boolean
     */
    boolean isStable(ICluster cluster);

    /**
     * if true the cluster is stable and will not allow removal
     *
     * @return Boolean
     */
    boolean isSemiStable(ICluster cluster);
}

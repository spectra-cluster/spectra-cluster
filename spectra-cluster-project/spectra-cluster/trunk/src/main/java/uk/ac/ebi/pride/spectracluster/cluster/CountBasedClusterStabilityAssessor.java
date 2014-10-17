package uk.ac.ebi.pride.spectracluster.cluster;

import uk.ac.ebi.pride.spectracluster.util.*;

/**
 * Cluster stability assessor based on the spectra count
 *
 * @author Steve Lewis
 * @author Rui Wang
 * @version $Id$
 */
public class CountBasedClusterStabilityAssessor implements IClusterStabilityAssessor {


    private final int stableClusterSize;
    private final int semiStableClusterSize;

    public CountBasedClusterStabilityAssessor() {
        this.stableClusterSize = StableClusterUtilities.getStableClusterSize();
        this.semiStableClusterSize = StableClusterUtilities.getSemiStableClusterSize();
    }

    public CountBasedClusterStabilityAssessor(int stableClusterSize, int semiStableClusterSize) {
        this.stableClusterSize = stableClusterSize;
        this.semiStableClusterSize = semiStableClusterSize;
    }

    @Override
    public boolean isStable(ICluster cluster) {
        int count = cluster.getClusteredSpectraCount();
        if (count == 1)
            return false; // Duh but saves other tests
        if (count >= stableClusterSize)
            return true;
          return false;
    }

    @Override
    public boolean isSemiStable(ICluster cluster) {
        int count = cluster.getClusteredSpectraCount();
        if (count == 1)
            return false; // Duh but saves other tests
        if (count >= semiStableClusterSize)
            return true;
        return false;
    }
}

package uk.ac.ebi.pride.spectracluster.cdf;

import uk.ac.ebi.pride.spectracluster.cluster.ICluster;

/**
 * This interface describes classes that derive the
 * number of comparisons to use for the threshold
 * calculation.
 *
 * Created by jg on 13.10.17.
 */
public interface INumberOfComparisonAssessor {
    /**
     * Returns the number of comparisons to use for the defined
     * cluster.
     * @param clusterToCompare The cluster to compare with the existing clusters.
     * @param nCurrentClusters The number of clusters this cluster will be compared to.
     * @return The number of comparisons to use.
     */
    int getNumberOfComparisons(ICluster clusterToCompare, int nCurrentClusters);
}

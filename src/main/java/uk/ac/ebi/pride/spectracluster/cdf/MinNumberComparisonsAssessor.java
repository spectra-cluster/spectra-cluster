package uk.ac.ebi.pride.spectracluster.cdf;

import uk.ac.ebi.pride.spectracluster.cluster.ICluster;

/**
 * Derives the number of comparisons to use based on the number
 * of existing clusters or the minimum number of comparisons to
 * use set in the constructor.
 *
 * Created by jg on 13.10.17.
 */
public class MinNumberComparisonsAssessor implements INumberOfComparisonAssessor {
    private final int minNumberComparisons;

    /**
     * Create an instante of the MinNumberComparisonsAssessor.
     * @param minNumberComparisons The minimum number of comparisons to use.
     */
    public MinNumberComparisonsAssessor(int minNumberComparisons) {
        this.minNumberComparisons = minNumberComparisons;
    }

    @Override
    public int getNumberOfComparisons(ICluster clusterToCompare, int nCurrentClusters) {
        if (nCurrentClusters < minNumberComparisons) {
            return minNumberComparisons;
        } else {
            return nCurrentClusters;
        }
    }

    public int getMinNumberComparisons() {
        return minNumberComparisons;
    }
}

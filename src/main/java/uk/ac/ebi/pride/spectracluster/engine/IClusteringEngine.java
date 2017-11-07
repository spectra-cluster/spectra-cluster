package uk.ac.ebi.pride.spectracluster.engine;

import uk.ac.ebi.pride.spectracluster.cluster.ICluster;
import uk.ac.ebi.pride.spectracluster.similarity.ISimilarityChecker;

import java.io.Serializable;
import java.util.Collection;

/**
 * This object does the clusters
 *
 * @author Johannes Griss
 * @author Steve Lewis
 * @author Rui Wang
 */
public interface IClusteringEngine extends Serializable{

    /**
     * Get the similarity check used
     *
     * @return an instance of similarity checker
     */
    ISimilarityChecker getSimilarityChecker();


    /**
     * Get the similarity threshold used
     *
     * @return similarity threshold
     */
    double getSimilarityThreshold();

    /**
     * Get clustered clusters sorted by MZ is useful
     *
     * @return !null list this will be sorted by mz a include clusters of all sizes
     */
    Collection<ICluster> getClusters();

    /**
     * add some clusters
     */
    void addClusters(ICluster... cluster);

    /**
     * clusters are merged in the internal collection
     *
     * @return true is  anything happened
     */
    boolean processClusters();

    /**
     * total number of clusters including queued clustersToAdd
     *
     * @return Integer representing the number of clusters
     */
    int size();
}

package uk.ac.ebi.pride.spectracluster.engine;

import uk.ac.ebi.pride.spectracluster.cluster.ICluster;

import java.util.Collection;

/**
 * uk.ac.ebi.pride.spectracluster.engine.IClusteringEngine
 * This object does the clusters
 *
 * @author Steve Lewis
 * @author Rui Wang
 */
public interface IIncrementalClusteringEngine extends IClusteringEngine {


    /**
     * add one cluster and return any clusters which are too far in mz from further consideration
     * NOTE clusters MUST be added in ascending MZ order
     *
     * @param added !null cluster to add
     * @return !null list of clusters not far enough away they will no longer change
     */
    Collection<ICluster> addClusterIncremental(ICluster added);


}

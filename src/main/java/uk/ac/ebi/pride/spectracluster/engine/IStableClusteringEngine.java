package uk.ac.ebi.pride.spectracluster.engine;

import uk.ac.ebi.pride.spectracluster.cluster.ICluster;

import java.util.Collection;

/**
 * @author Steve Lewis
 * @author Rui Wang
 * @version $Id$
 *          <p/>
 *          todo: development
 */
@Deprecated
public interface IStableClusteringEngine extends IIncrementalClusteringEngine {

    public void addUnstableCluster(ICluster unstableCluster);

    public void processStableCluster(ICluster stableCluster);

    public Collection<ICluster> getClusters();
}

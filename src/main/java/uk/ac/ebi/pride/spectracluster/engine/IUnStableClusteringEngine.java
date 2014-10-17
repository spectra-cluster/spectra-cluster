package uk.ac.ebi.pride.spectracluster.engine;

import uk.ac.ebi.pride.spectracluster.cluster.ICluster;

/**
 * @author Steve Lewis
 * @author Rui Wang
 * @version $Id$
 *          <p/>
 *          todo: development
 */
@Deprecated
public interface IUnStableClusteringEngine {

    public void addStableCluster(ICluster unstableCluster);

    /**
     * try to move spectra to the stable cluster
     *
     * @param unstableCluster
     * @return true if changed
     */
    public boolean processUnStableCluster(ICluster unstableCluster);
}

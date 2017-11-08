package uk.ac.ebi.pride.spectracluster.util;

import uk.ac.ebi.pride.spectracluster.cluster.ICluster;

import java.io.Serializable;

/**
 * uk.ac.ebi.pride.spectracluster.util.ClusterCreateListener
 * make a listener for cluster reads
 * User: Steve
 * Date: 9/23/13
 */
public interface ClusterCreateListener extends Serializable{

    /**
     * initialize reading - if reading happens once - say from
     * one file all this may happen in the constructor
     */
    void onClusterStarted(Object... otherData);

    /**
     * do something when a cluster is created or read
     *
     * @param cluster
     */
    void onClusterCreate(ICluster cluster, Object... otherData);

    /**
     * do something when a cluster when the last cluster is read -
     * this may be after a file read is finished
     */
    void onClusterCreateFinished(Object... otherData);


}

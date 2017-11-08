package uk.ac.ebi.pride.spectracluster.util;

import java.io.Serializable;

/**
 * uk.ac.ebi.pride.spectracluster.util.IAlgorithm
 * many interfaces are algorithms suhc as clusteringenging, concensusspectrum ...
 * This marks them as algorithms and forced them to have names and versions
 * User: Steve
 * Date: 7/17/13
 */
public interface IAlgorithm extends Serializable {

    /**
     * return a name which should not change
     *
     * @return !null name
     */
    public String getName();

    /**
     * return a version number - this may be updated over time
     *
     * @return !null version
     */
    public String getCurrentVersion();

}

package uk.ac.ebi.pride.spectracluster.io;

import uk.ac.ebi.pride.spectracluster.cluster.ICluster;

/**
 * uk.ac.ebi.pride.spectracluster.cluster.IClusterWriter
 * Abstract the concept of appending a cluster to an appendable.
 * This may include filter
 * User: Steve
 * Date: 9/25/13
 */
public interface IClusterAppender {
    /**
     * @param out       !null open appendable
     * @param data      !null cluster
     * @param otherData any other data - implementation specific and usually blank
     * @return true if anything was appended otherwise false
     */
    void appendCluster(Appendable out, ICluster data, Object... otherData) throws AppenderException;

    /**
     * add whatever happens at the start
     *
     * @param out       !null open appendale
     * @param otherData any other data - implementation specific and usually blank
     * @return true if anything was appended otherwise false
     */
    void appendStart(Appendable out, Object... otherData) throws AppenderException;

    /**
     * add whatever happens at the end
     *
     * @param out       !null open appendale
     * @param otherData any other data - implementation specific and usually blank
     * @return true if anything was appended otherwise false
     */
    void appendEnd(Appendable out, Object... otherData) throws AppenderException;
}

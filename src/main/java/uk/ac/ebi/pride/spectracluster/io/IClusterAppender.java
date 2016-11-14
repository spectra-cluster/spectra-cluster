package uk.ac.ebi.pride.spectracluster.io;

import uk.ac.ebi.pride.spectracluster.cluster.ICluster;

/**
 * Abstract the concept of appending a cluster to an appendable.
 *
 * @author Johannes Griss
 */
public interface IClusterAppender {
    /**
     * @param out       !null open appendable
     * @param data      !null cluster
     * @param otherData any other data - implementation specific and usually blank
     */
    void appendCluster(Appendable out, ICluster data, Object... otherData) throws AppenderException;

    /**
     * add whatever happens at the start
     *
     * @param out       !null open appendale
     * @param otherData any other data - implementation specific and usually blank
     */
    void appendStart(Appendable out, Object... otherData) throws AppenderException;

    /**
     * add whatever happens at the end
     *
     * @param out       !null open appendale
     * @param otherData any other data - implementation specific and usually blank
     */
    void appendEnd(Appendable out, Object... otherData) throws AppenderException;
}

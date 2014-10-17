package uk.ac.ebi.pride.spectracluster.io;

import uk.ac.ebi.pride.spectracluster.cluster.ICluster;

/**
 * uk.ac.ebi.pride.spectracluster.cluster.CGFClusterAppender
 * User: Steve
 * Date: 9/25/13
 */

public class MSFClusterAppender implements IClusterAppender {

    public static MSFClusterAppender INSTANCE = new MSFClusterAppender(MSFSpectrumAppender.INSTANCE);

    private final MSFSpectrumAppender spectrumAppender;

    public MSFClusterAppender(MSFSpectrumAppender spectrumAppender) {
        this.spectrumAppender = spectrumAppender;
    }

    /**
     * @param out       !null open appendale
     * @param data      !null cluster
     * @param otherData any other data - implementation specific and usually blank
     * @return true if anything was appended otherwise false
     */
    @Override
    public void appendCluster(final Appendable out, final ICluster data, final Object... otherData) {
        spectrumAppender.appendSpectrum(out, data.getConsensusSpectrum());
    }

    /**
     * add whatever happens at the start
     *
     * @param out       !null open appendale
     * @param otherData any other data - implementation specific and usually blank
     * @return true if anything was appended otherwise false
     */
    @Override
    public void appendStart(final Appendable out, final Object... otherData) {
    }

    /**
     * add whatever happens at the end
     *
     * @param out       !null open appendale
     * @param otherData any other data - implementation specific and usually blank
     * @return true if anything was appended otherwise false
     */
    @Override
    public void appendEnd(final Appendable out, final Object... otherData) {
    }
}

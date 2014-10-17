package uk.ac.ebi.pride.spectracluster.io;

import uk.ac.ebi.pride.spectracluster.cluster.ICluster;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;

import java.io.IOException;
import java.util.List;

/**
 * uk.ac.ebi.pride.spectracluster.cluster.CGFClusterAppender
 * User: Steve
 * Date: 9/25/13
 */
public class CGFClusterAppender implements IClusterAppender {


    public static CGFClusterAppender INSTANCE = new CGFClusterAppender(MGFSpectrumAppender.INSTANCE);



    private final MGFSpectrumAppender spectrumAppender;

    private CGFClusterAppender(MGFSpectrumAppender spectrumAppender) {
        this.spectrumAppender = spectrumAppender;
    }

    /**
     * @param out       !null open appendale
     * @param cluster   !null cluster
     * @param otherData any other cluster - implementation specific and usually blank
     */
    @Override
    public void appendCluster(final Appendable out, final ICluster cluster, final Object... otherData) {
        try {
            out.append("BEGIN CLUSTER");
            out.append(" Id=").append(cluster.getId());
            out.append(" Charge=").append(String.valueOf(cluster.getPrecursorCharge()));

            out.append("\n");

            appendSpectra(out, cluster);

            out.append("END CLUSTER");
            out.append("\n");
        } catch (IOException e) {
            throw new AppenderException(e);
        }
    }

    private void appendSpectra(final Appendable out, final ICluster cluster) {
        List<ISpectrum> clusteredSpectra = cluster.getClusteredSpectra();
        for (ISpectrum cs : clusteredSpectra) {
            spectrumAppender.appendSpectrum(out, cs);  // single spectrum become mgfs

        }
    }


    /**
     * add whatever happens at the start
     *
     * @param out       !null open appendale
     * @param otherData any other data - implementation specific and usually blank
     */
    @Override
    public void appendStart(final Appendable out, final Object... otherData) {
    }

    /**
     * add whatever happens at the end
     *
     * @param out       !null open appendale
     * @param otherData any other data - implementation specific and usually blank
     */
    @Override
    public void appendEnd(final Appendable out, final Object... otherData) {
    }
}

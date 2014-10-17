package uk.ac.ebi.pride.spectracluster.io;

import uk.ac.ebi.pride.spectracluster.cluster.ICluster;

import java.io.IOException;

/**
 * Append cluster into TSV file
 *
 * @author Rui Wang
 * @version $Id$
 */
public class TSVClusterAppender implements IClusterAppender {

    @Override
    public void appendCluster(Appendable out, ICluster cluster, Object... otherData) {
        try {
            out.append(cluster.getId());
            out.append("\t");

            String mz = String.format("%f8.2", cluster.getPrecursorMz());
            out.append(mz);
            out.append("\t");

            out.append(Integer.toString(cluster.getPrecursorCharge()));
            out.append("\t");
        } catch (IOException e) {
            throw new AppenderException(e);
        }
    }

    @Override
    public void appendStart(Appendable out, Object... otherData) {
    }

    @Override
    public void appendEnd(Appendable out, Object... otherData) {
    }
}

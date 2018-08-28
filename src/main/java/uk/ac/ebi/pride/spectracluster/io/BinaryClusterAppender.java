package uk.ac.ebi.pride.spectracluster.io;

import uk.ac.ebi.pride.spectracluster.cluster.ICluster;
import uk.ac.ebi.pride.spectracluster.consensus.IConsensusSpectrumBuilder;
import uk.ac.ebi.pride.spectracluster.spectrum.IPeak;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.util.ComparisonMatch;

import java.io.IOException;
import java.io.ObjectOutputStream;
import java.util.List;
import java.util.Properties;

/**
 * uk.ac.ebi.pride.spectracluster.cluster.CGFClusterAppender
 * User: Steve
 * Date: 9/25/13
 */
public class BinaryClusterAppender {


    public static BinaryClusterAppender INSTANCE = new BinaryClusterAppender();

    private BinaryClusterAppender() {

    }

    /**
     * @param out       !null open ObjectOutputStream
     * @param cluster   !null cluster
     */
    public void appendCluster(final ObjectOutputStream out, final ICluster cluster) {
        try {
            // first save the classname
            out.writeObject(cluster.getClass().getCanonicalName());

            // standard fields
            out.writeObject(cluster.getId());
            out.writeInt(cluster.getPrecursorCharge());
            out.writeFloat(cluster.getPrecursorMz());
            out.writeBoolean(cluster.storesPeakLists());

            appendComparisonMatches(out, cluster);

            if (!cluster.storesPeakLists())
                appendConsensusSpectrumBuilder(out, cluster.getConsensusSpectrumBuilder());

            appendSpectra(out, cluster);
        } catch (IOException e) {
            throw new AppenderException(e);
        }
    }

    private void appendConsensusSpectrumBuilder(ObjectOutputStream out, IConsensusSpectrumBuilder consensusSpectrumBuilder) throws IOException {
        // always save the class first
        out.writeObject(consensusSpectrumBuilder.getClass().getCanonicalName());

        // standard fields
        out.writeObject(consensusSpectrumBuilder.getConsensusSpectrum().getId());
        out.writeInt(consensusSpectrumBuilder.getSpectraCount());
        out.writeInt(consensusSpectrumBuilder.getSumCharge());
        out.writeDouble(consensusSpectrumBuilder.getSumPrecursorIntensity());
        out.writeDouble(consensusSpectrumBuilder.getSumPrecursorMz());

        // write the raw peak list
        appendPeaklist(out, consensusSpectrumBuilder.getRawConsensusPeaks());
    }

    private void appendPeaklist(ObjectOutputStream out, List<IPeak> peaklist) throws IOException {
        out.writeInt(peaklist.size());

        for (IPeak peak : peaklist) {
            out.writeFloat(peak.getMz());
            out.writeFloat(peak.getIntensity());
            out.writeInt(peak.getCount());
        }
    }

    private void appendComparisonMatches(ObjectOutputStream out, ICluster cluster) throws IOException {
        out.writeInt(cluster.getComparisonMatches().size());

        // ignore empty matches
        if (cluster.getComparisonMatches().isEmpty())
            return;

        // write the comparison matches
        for (ComparisonMatch c : cluster.getComparisonMatches()) {
            out.writeFloat(c.getSimilarity());
            out.writeObject(c.getSpectrumId());
        }
    }

    private void appendSpectra(ObjectOutputStream out, ICluster cluster) throws IOException {
        List<ISpectrum> clusteredSpectra = cluster.getClusteredSpectra();
        out.writeInt(clusteredSpectra.size());

        for (ISpectrum cs : clusteredSpectra) {
            // default properties
            out.writeObject(cs.getId());
            out.writeInt(cs.getPrecursorCharge());
            out.writeFloat(cs.getPrecursorMz());

            // additional properties
            Properties properties = cs.getProperties();
            out.writeObject(properties);

            // peak list
            appendPeaklist(out, cs.getPeaks());
        }
    }

    public void appendEnd(ObjectOutputStream out) throws IOException {
        out.writeObject("END");
    }
}

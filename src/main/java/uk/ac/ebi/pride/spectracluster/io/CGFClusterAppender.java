package uk.ac.ebi.pride.spectracluster.io;

import uk.ac.ebi.pride.spectracluster.cluster.ICluster;
import uk.ac.ebi.pride.spectracluster.consensus.IConsensusSpectrumBuilder;
import uk.ac.ebi.pride.spectracluster.spectrum.IPeak;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.util.ComparisonMatch;

import java.io.IOException;
import java.util.List;
import java.util.Properties;

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
            out.append(" ContainsPeaklist=").append(String.valueOf(cluster.storesPeakLists()));

            out.append("\n");

            appendProperties(out, cluster);
            appendComparisonMatches(out, cluster);

            if (!cluster.storesPeakLists())
                appendConsensusSpectrumBuilder(out, cluster.getConsensusSpectrumBuilder());

            appendSpectra(out, cluster);

            out.append("END CLUSTER");
            out.append("\n");
        } catch (IOException e) {
            throw new AppenderException(e);
        }
    }

    private void appendProperties(Appendable out, ICluster cluster) throws IOException {
        out.append("Properties=");
        Properties properties = cluster.getProperties();

        boolean firstProperty = true;

        for (String name : properties.stringPropertyNames()) {
            if (!firstProperty)
                out.append("#");
            out.append(name + "=" + properties.getProperty(name));
        }

        out.append("\n");
    }

    private void appendConsensusSpectrumBuilder(Appendable out, IConsensusSpectrumBuilder consensusSpectrumBuilder) throws IOException {
        StringBuilder consensusSpectrumString = new StringBuilder();

        consensusSpectrumString.append("BEGIN CONSENSUS");
        consensusSpectrumString.append(" id=").append(consensusSpectrumBuilder.getConsensusSpectrum().getId());
        consensusSpectrumString.append(" class=").append(consensusSpectrumBuilder.getClass().getCanonicalName());
        consensusSpectrumString.append(" nSpec=").append(consensusSpectrumBuilder.getSpectraCount());
        consensusSpectrumString.append(" SumCharge=").append(consensusSpectrumBuilder.getSumCharge());
        consensusSpectrumString.append(" SumIntens=").append(consensusSpectrumBuilder.getSumPrecursorIntensity());
        consensusSpectrumString.append(" SumMz=").append(consensusSpectrumBuilder.getSumPrecursorMz());
        consensusSpectrumString.append("\n");

        for (IPeak peak : consensusSpectrumBuilder.getRawConsensusPeaks()) {
            String line = String.format("%10.3f", peak.getMz()).trim() + "\t" +
                    String.format("%10.3f", peak.getIntensity()).trim() + "\t" + peak.getCount();
            consensusSpectrumString.append(line).append("\n");
        }

        consensusSpectrumString.append("END CONSENSUS\n");

        out.append(consensusSpectrumString.toString());
    }

    private void appendComparisonMatches(Appendable out, ICluster cluster) throws IOException {
        // ignore empty matches
        if (cluster.getComparisonMatches().isEmpty())
            return;

        StringBuilder comparisonMatchesString = new StringBuilder("ComparisonMatches=");

        for (ComparisonMatch c : cluster.getComparisonMatches()) {
            if (comparisonMatchesString.length() > 19)
                comparisonMatchesString.append("#");

            comparisonMatchesString.append(c.getSimilarity()).append(":").append(c.getSpectrumId());
        }

        comparisonMatchesString.append("\n");

        out.append(comparisonMatchesString.toString());
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

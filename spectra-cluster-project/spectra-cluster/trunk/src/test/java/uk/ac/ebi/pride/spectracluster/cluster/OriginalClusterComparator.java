package uk.ac.ebi.pride.spectracluster.cluster;

import uk.ac.ebi.pride.spectracluster.quality.IQualityScorer;

import java.util.Comparator;

/**
 * @author Rui Wang
 * @version $Id$
 */
public class OriginalClusterComparator implements Comparator<ICluster> {
    private static final double PRECURSOR_RANGE = 2;

    private IQualityScorer qualityScorer;

    public OriginalClusterComparator(IQualityScorer qualityScorer) {
        this.qualityScorer = qualityScorer;
    }

    @Override
    public int compare(ICluster cluster1, ICluster cluster2) {
        // first check whether the precursor m/z is different
        if (cluster1.getPrecursorMz() - PRECURSOR_RANGE > cluster2.getPrecursorMz())
            return -1;
        if (cluster1.getPrecursorMz() + PRECURSOR_RANGE < cluster2.getPrecursorMz())
            return 1;

        // as the m/z ranges are the same check the quality
        double quality1 = qualityScorer.calculateQualityScore(cluster1.getConsensusSpectrum());
        double quality2 = qualityScorer.calculateQualityScore(cluster2.getConsensusSpectrum());

        if (quality1 > quality2)
            return -1;
        if (quality1 < quality2)
            return 1;

        return 0;
    }
}

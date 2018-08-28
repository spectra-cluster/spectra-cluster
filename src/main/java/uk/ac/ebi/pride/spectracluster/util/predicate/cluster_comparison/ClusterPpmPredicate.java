package uk.ac.ebi.pride.spectracluster.util.predicate.cluster_comparison;

import uk.ac.ebi.pride.spectracluster.cluster.ICluster;
import uk.ac.ebi.pride.spectracluster.util.predicate.IComparisonPredicate;

/**
 * This predicate is used to test whether two clusters
 * are within a certain precursor ppm window.
 *
 * Created by jg on 24.06.17.
 */
public class ClusterPpmPredicate implements IComparisonPredicate<ICluster> {
    private final float mzFraction;

    /**
     * Constructs a new ClusterPpmPredicate that tests whether the two
     * passed clusters' precursor m/z are within a certain ppm tolerance
     *
     * @param ppmThreshold The threshold in ppm.
     */
    public ClusterPpmPredicate(float ppmThreshold) {
        this.mzFraction = ppmThreshold / 1000000;
    }

    @Override
    public boolean apply(ICluster o1, ICluster o2) {
        // base the tolerance on the higher m/z value
        float maxMz = o1.getPrecursorMz() > o2.getPrecursorCharge() ? o1.getPrecursorMz() : o2.getPrecursorMz();

        // calculate the tolerance in m/z
        float maxDifference = maxMz * this.mzFraction;

        // get the difference
        float difference = Math.abs(o1.getPrecursorMz() - o2.getPrecursorMz());

        return difference <= maxDifference;
    }
}

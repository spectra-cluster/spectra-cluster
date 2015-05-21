package uk.ac.ebi.pride.spectracluster.util.predicate.cluster_comparison;

import uk.ac.ebi.pride.spectracluster.cluster.ICluster;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.util.predicate.IComparisonPredicate;
import uk.ac.ebi.pride.spectracluster.util.predicate.spectrum_comparison.ShareMajorPeaksPredicate;

/**
 * Created by jg on 20.05.15.
 */
public class ClusterShareMajorPeakPredicate implements IComparisonPredicate<ICluster> {
    private final IComparisonPredicate<ISpectrum> majorPeakPredicate;

    public ClusterShareMajorPeakPredicate() {
        this.majorPeakPredicate = new ShareMajorPeaksPredicate();
    }

    public ClusterShareMajorPeakPredicate(int nMajorPeaks) {
        this.majorPeakPredicate = new ShareMajorPeaksPredicate(nMajorPeaks);
    }

    @Override
    public boolean apply(ICluster o1, ICluster o2) {
        return majorPeakPredicate.apply(o1.getConsensusSpectrum(), o2.getConsensusSpectrum());
    }
}

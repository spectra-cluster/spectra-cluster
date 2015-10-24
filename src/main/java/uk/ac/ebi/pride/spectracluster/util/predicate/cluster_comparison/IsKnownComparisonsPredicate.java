package uk.ac.ebi.pride.spectracluster.util.predicate.cluster_comparison;

import uk.ac.ebi.pride.spectracluster.cluster.ICluster;
import uk.ac.ebi.pride.spectracluster.util.predicate.IComparisonPredicate;

/**
 * Created by jg on 20.05.15.
 */
public class IsKnownComparisonsPredicate implements IComparisonPredicate<ICluster> {
    @Override
    public boolean apply(ICluster o1, ICluster o2) {
        // make sure both have comparison matches
        if (o1.getComparisonMatches().size() == 0 || o2.getComparisonMatches().size() == 0)
            return false;

        // check if they are known
        if (o1.isKnownComparisonMatch(o2.getId()) || o2.isKnownComparisonMatch(o1.getId()))
            return true;

        // not known
        return false;
    }
}

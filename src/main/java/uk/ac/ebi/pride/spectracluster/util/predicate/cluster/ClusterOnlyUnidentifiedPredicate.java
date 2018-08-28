package uk.ac.ebi.pride.spectracluster.util.predicate.cluster;

import uk.ac.ebi.pride.spectracluster.cluster.ICluster;
import uk.ac.ebi.pride.spectracluster.spectrum.KnownProperties;
import uk.ac.ebi.pride.spectracluster.util.predicate.IPredicate;

/**
 * Only applies to cluster that contain at only unidentified spectrum.
 *
 * Created by jg on 14.01.18.
 */
public class ClusterOnlyUnidentifiedPredicate implements IPredicate<ICluster> {
    @Override
    public boolean apply(ICluster o) {
        return o.getClusteredSpectra()
                .stream()
                .allMatch(spectrum -> spectrum.getProperty(KnownProperties.IDENTIFIED_PEPTIDE_KEY) == null);
    }
}

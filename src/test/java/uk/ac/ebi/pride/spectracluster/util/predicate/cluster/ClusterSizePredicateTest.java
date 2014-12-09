package uk.ac.ebi.pride.spectracluster.util.predicate.cluster;

import org.junit.Test;
import uk.ac.ebi.pride.spectracluster.cluster.ICluster;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import static org.mockito.Mockito.mock;
import static org.mockito.Mockito.when;

public class ClusterSizePredicateTest {

    @Test
    public void testApply() throws Exception {
        ClusterSizePredicate clusterSizePredicate = new ClusterSizePredicate();

        ICluster cluster = mock(ICluster.class);
        when(cluster.getClusteredSpectraCount()).thenReturn(10);
        assertTrue(clusterSizePredicate.apply(cluster));


        when(cluster.getClusteredSpectraCount()).thenReturn(9);
        assertFalse(clusterSizePredicate.apply(cluster));
    }
}
package uk.ac.ebi.pride.spectracluster.util.predicate.cluster_comparison;

import junit.framework.Assert;
import org.junit.Before;
import org.junit.Test;
import uk.ac.ebi.pride.spectracluster.cluster.ICluster;
import uk.ac.ebi.pride.spectracluster.io.ParserUtilities;
import uk.ac.ebi.pride.spectracluster.util.predicate.IComparisonPredicate;

import java.io.File;
import java.util.List;

/**
 * Created by jg on 24.06.17.
 */
public class ClusterPpmPredicateTest {
    protected File testFile;

    @Before
    public void setUp() throws Exception {
        testFile = new File(ClusterPpmPredicateTest.class.getClassLoader().getResource("kuester_test.mgf").toURI());
    }

    @Test
    public void testPredicate() throws Exception {
        List<ICluster> clusters = ParserUtilities.readMGFClusters(testFile);

        IComparisonPredicate<ICluster> predicate = new ClusterPpmPredicate(20);

        Assert.assertTrue(predicate.apply(clusters.get(0), clusters.get(1)));
        Assert.assertFalse(predicate.apply(clusters.get(0), clusters.get(clusters.size() - 1)));
    }
}

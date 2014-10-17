package uk.ac.ebi.pride.spectracluster.cluster;

import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;
import uk.ac.ebi.pride.spectracluster.engine.*;
import uk.ac.ebi.pride.spectracluster.util.ClusteringTestUtilities;
import uk.ac.ebi.pride.spectracluster.util.Defaults;
import uk.ac.ebi.pride.spectracluster.util.StableClusterUtilities;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

/**
 * @author Rui Wang
 * @version $Id$
 */
public class StableClusteringEngineTest {

    public static final int STABLE_CLUSTER_SIZE = 4;

    private List<ICluster> unstableClusters;
    private List<ICluster> stableClusters;
    private StableClusteringEngine clusteringEngine;

    @Before
    public void setUp() throws Exception {
        StableClusterUtilities.setStableClusterSize(STABLE_CLUSTER_SIZE);   // drop cluster size foe a small sample
        IClusteringEngine ce = EngineFactories.DEFAULT_CLUSTERING_ENGINE_FACTORY.buildInstance();
        List<ICluster> originalSpectralClusters = new ArrayList<ICluster>(ClusteringTestUtilities.readSpectraClustersFromResource());

        for (ICluster sc : originalSpectralClusters) {
            ce.addClusters(sc);
        }
        originalSpectralClusters = new ArrayList<ICluster>(ce.getClusters());

        clusteringEngine = new StableClusteringEngine(Defaults.getDefaultSimilarityChecker(), Defaults.getSimilarityThreshold());
        unstableClusters = new ArrayList<ICluster>();
        stableClusters = new ArrayList<ICluster>();

        final CountBasedClusterStabilityAssessor clusterStabilityAssessor = new CountBasedClusterStabilityAssessor(StableClusterUtilities.getStableClusterSize(), StableClusterUtilities.getSemiStableClusterSize());
        for (ICluster originalSpectralCluster : originalSpectralClusters) {
            if (clusterStabilityAssessor.isStable(originalSpectralCluster)) {
                stableClusters.add(originalSpectralCluster);
            } else {
                unstableClusters.add(originalSpectralCluster);
            }
        }
    }


    public static final int NUMBER_MERGED_CLUSTERS = 0;

    @Test
    public void testEngine() throws Exception {
        Assert.assertTrue(!unstableClusters.isEmpty());
        Assert.assertTrue(!stableClusters.isEmpty());

        for (ICluster unstableCluster : unstableClusters) {
            clusteringEngine.addUnstableCluster(unstableCluster);
        }

        for (ICluster stableCluster : stableClusters) {
            int clusteredSpectraCount = stableCluster.getClusteredSpectraCount();
            int numberOfUnstableSpectra = clusteringEngine.getNumberOfUnstableSpectra();
            clusteringEngine.processStableCluster(stableCluster);
            int processedNumberOfUnstableSpectra = clusteringEngine.getNumberOfUnstableSpectra();
            int processedClusteredSpectrumCount = stableCluster.getClusteredSpectraCount();
//            Assert.assertEquals(clusteredSpectraCount,processedClusteredSpectrumCount);

            if (clusteredSpectraCount != processedClusteredSpectrumCount) {
                int actual = numberOfUnstableSpectra - processedNumberOfUnstableSpectra;
                int expected = processedClusteredSpectrumCount - clusteredSpectraCount;
                Assert.assertEquals(expected, actual);
            }
        }

        Collection<ICluster> unstableClustersAfterProcess = clusteringEngine.getClusters();
        // we did get some merging and I think that is good
        int actual = unstableClusters.size() - NUMBER_MERGED_CLUSTERS;
        Assert.assertEquals(unstableClustersAfterProcess.size(), actual);
    }
}

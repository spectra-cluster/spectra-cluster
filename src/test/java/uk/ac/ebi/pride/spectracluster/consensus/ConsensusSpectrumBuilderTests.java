package uk.ac.ebi.pride.spectracluster.consensus;

import org.junit.*;
import uk.ac.ebi.pride.spectracluster.cluster.ICluster;
import uk.ac.ebi.pride.spectracluster.filter.*;
import uk.ac.ebi.pride.spectracluster.spectrum.*;
import uk.ac.ebi.pride.spectracluster.util.*;

import java.util.List;

/**
 * uk.ac.ebi.pride.spectracluster.consensus.ConsensusSpectrumBuilderTests
 *
 * @author Steve Lewis
 *         This tests the difference between  IConsensusSpectrumBuilder implementations
 */
public class ConsensusSpectrumBuilderTests {


    @Test
    public void testConsensusSpectrum() throws Exception {

        // do nothing to filter
        Defaults.setDefaultPeakFilter(IPeakFilter.NULL_FILTER);


        IConsensusSpectrumBuilder factory1 = ConsensusSpectrum.buildFactory(IPeakFilter.NULL_FILTER).getConsensusSpectrumBuilder();
        IConsensusSpectrumBuilder factory2 = new JohannesConsensusSpectrum();


        List<ICluster> clusters = ClusteringTestUtilities.readSpectraClustersFromResource();

        List<ISpectrum> fromFactory1 = ClusteringTestUtilities.buildConsessusSpectra(clusters, factory1);

        List<ISpectrum> fromFactory2 = ClusteringTestUtilities.buildConsessusSpectra(clusters, factory2);


        Assert.assertEquals(fromFactory1.size(), fromFactory2.size());
        for (int i = 0; i < fromFactory1.size(); i++) {
            final ISpectrum oldSpec = fromFactory2.get(i);
            final ISpectrum newSpec = fromFactory1.get(i);
            boolean equivalent = ClusteringTestUtilities.areSpectraVeryClose(newSpec, oldSpec);
            final double v = Defaults.getDefaultSimilarityChecker().assessSimilarity(oldSpec, newSpec);
            if (!equivalent) {

                equivalent = ClusteringTestUtilities.areSpectraVeryClose(newSpec, oldSpec);// break here do debug failure
                Assert.assertTrue(equivalent);
            }

        }

    }


}

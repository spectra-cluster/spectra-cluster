package uk.ac.ebi.pride.spectracluster.util;

import junit.framework.Assert;
import org.junit.Test;
import uk.ac.ebi.pride.spectracluster.consensus.ConsensusSpectrum;
import uk.ac.ebi.pride.spectracluster.consensus.GreedyConsensusSpectrum;
import uk.ac.ebi.pride.spectracluster.consensus.IConsensusSpectrumBuilder;
import uk.ac.ebi.pride.spectracluster.engine.IClusteringEngine;
import uk.ac.ebi.pride.spectracluster.similarity.ISimilarityChecker;

/**
 * Tests whether setting of default parameters propagates through
 * other classes
 *
 * Created by jg on 07.05.17.
 */
public class DefaultsTest {
    @Test
    public void testFragmentIonTolerance() throws Exception {
        float targetTolerance = 0.123F;
        Defaults.setFragmentIonTolerance(targetTolerance);

        IConsensusSpectrumBuilder consensusSpectrumBuilder = GreedyConsensusSpectrum.buildFactory().getConsensusSpectrumBuilder();
        Assert.assertEquals(targetTolerance, consensusSpectrumBuilder.getFragmentIonTolerance());

        consensusSpectrumBuilder = ConsensusSpectrum.buildFactory().getConsensusSpectrumBuilder();
        Assert.assertEquals(targetTolerance, consensusSpectrumBuilder.getFragmentIonTolerance());

        // test different other default classes that use this tolerance
        ISimilarityChecker similarityChecker = Defaults.getDefaultSimilarityChecker();
        Assert.assertEquals(targetTolerance, similarityChecker.getFragmentIonTolerance());

        IClusteringEngine clusteringEngine = Defaults.getDefaultClusteringEngine();
        Assert.assertEquals(targetTolerance, clusteringEngine.getSimilarityChecker().getFragmentIonTolerance());
    }
}

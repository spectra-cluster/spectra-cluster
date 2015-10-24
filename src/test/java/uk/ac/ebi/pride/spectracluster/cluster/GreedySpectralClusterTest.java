package uk.ac.ebi.pride.spectracluster.cluster;

import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;
import uk.ac.ebi.pride.spectracluster.consensus.ConsensusSpectrum;
import uk.ac.ebi.pride.spectracluster.consensus.IConsensusSpectrumBuilder;
import uk.ac.ebi.pride.spectracluster.io.ParserUtilities;
import uk.ac.ebi.pride.spectracluster.similarity.FrankEtAlDotProduct;
import uk.ac.ebi.pride.spectracluster.similarity.ISimilarityChecker;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.util.function.peak.NullPeakFunction;

import java.io.File;
import java.util.List;

/**
 * Created by jg on 11.05.15.
 */
public class GreedySpectralClusterTest {
    private ISpectrum[] testSpectra;

    @Before
    public void setUp() throws Exception {
        File testFile = new File(GreedySpectralClusterTest.class.getClassLoader().getResource("uk/ac/ebi/pride/spectracluster/engine/cluster_spec_4776.mgf").toURI());
        testSpectra = ParserUtilities.readMGFScans(testFile);
    }

    @Test
    public void testAddSpectra() throws Exception {
        GreedySpectralCluster spectralCluster = new GreedySpectralCluster("testId");

        spectralCluster.addSpectra(testSpectra);
        Assert.assertEquals(testSpectra.length, spectralCluster.getClusteredSpectraCount());
    }

    @Test
    public void testConsensusSpectrum() throws Exception {
        GreedySpectralCluster spectralCluster = new GreedySpectralCluster("testId");

        spectralCluster.addSpectra(testSpectra);
        ISpectrum greedyConsensusSpectrum = spectralCluster.getConsensusSpectrum();

        // build a "standard" consensus spectrum as point of reference
        IConsensusSpectrumBuilder referenceConsensusSpectrum = ConsensusSpectrum.buildFactory().getConsensusSpectrumBuilder();
        referenceConsensusSpectrum.addSpectra(testSpectra);
        ISpectrum referenceSpec = referenceConsensusSpectrum.getConsensusSpectrum();

        ISimilarityChecker similarityChecker = new FrankEtAlDotProduct(0.5F);
        similarityChecker.setPeakFiltering(true);
        double similarity = similarityChecker.assessSimilarity(greedyConsensusSpectrum, referenceSpec);

        // spectra must be nearly identical
        Assert.assertEquals(0.953, similarity, 0.001);
    }

    @Test
    public void testSaveComparisonMatches() {
        double[] similarities = {0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.4, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.5, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.6, 0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68, 0.69, 0.7, 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 0.8, 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1};

        GreedySpectralCluster spectralCluster = new GreedySpectralCluster("test");

        for (int i = 0; i < similarities.length; i++) {
            spectralCluster.saveComparisonResult(String.valueOf(i), (float) similarities[i]);
        }

        Assert.assertEquals(30, spectralCluster.getComparisonMatches().size());
        Assert.assertFalse(spectralCluster.isInBestComparisonResults("1"));
        Assert.assertFalse(spectralCluster.isInBestComparisonResults("60"));
        Assert.assertTrue(spectralCluster.isInBestComparisonResults("61"));
    }
}

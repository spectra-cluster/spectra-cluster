package uk.ac.ebi.pride.spectracluster.quality;

import junit.framework.Assert;
import org.junit.*;
import org.junit.Test;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.util.*;

import java.util.List;

/**
 * @author Rui Wang
 * @version $Id$
 */
public class SignalToNoiseQualityCheckerTests {
    private List<ISpectrum> peptideSpectrumMatches;
    private IQualityScorer originalQualityScorer;
    private IQualityScorer qualityScorer;


    @Before
    public void setUp() throws Exception {
        originalQualityScorer = new OriginalSignalToNoiseChecker();
        qualityScorer =  Defaults.getDefaultQualityScorer();
        peptideSpectrumMatches = ClusteringTestUtilities.readConsensusSpectralItems();
    }

    @Test
    public void testSignalToNoiseQualityChecker() throws Exception {
        for (ISpectrum peptideSpectrumMatch : peptideSpectrumMatches) {

            double originalScore = originalQualityScorer.calculateQualityScore(peptideSpectrumMatch);
            double score = qualityScorer.calculateQualityScore(peptideSpectrumMatch);
            Assert.assertEquals(score, originalScore, 0.1);
        }
    }
}

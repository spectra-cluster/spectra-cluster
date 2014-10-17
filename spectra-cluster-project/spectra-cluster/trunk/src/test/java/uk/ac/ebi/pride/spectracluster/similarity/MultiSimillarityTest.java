package uk.ac.ebi.pride.spectracluster.similarity;

import org.junit.Before;
import org.junit.Test;
import uk.ac.ebi.pride.spectracluster.io.ParserUtilities;
import uk.ac.ebi.pride.spectracluster.util.ConsensusSpectraItems;

import java.io.File;
import java.net.URL;


/**
 * @author Rui Wang
 * @version $Id$
 */
public class MultiSimillarityTest {
    private ConsensusSpectraItems[] consensusSpectraItems;

    @SuppressWarnings("UnusedDeclaration")
    public static final double SIMILAR_THRESHOLD = 1.5; // this is really big

    @Before
    public void setUp() throws Exception {
        // load a file contains a list of clusters
        // URL url = MultiSimillarityTest.class.getClassLoader().getResource ("uk/ac/ebi/pride/spectracluster/util/spectra_400.0_4.0.cgf");
        URL url = MultiSimillarityTest.class.getClassLoader().getResource("uk/ac/ebi/pride/spectracluster/io/spectra_400.0_4.0.cgf");
        if (url == null) {
            throw new IllegalStateException("no file for input found!");
        }
        File inputFile = new File(url.toURI());

        consensusSpectraItems = ParserUtilities.readClusters(inputFile);


    }

    /**
     * This test has been disabled since it is trying to use an non-optimized version
     * of the FrankEtAlDotProduct
     *
     * @throws Exception
     */
    @Test
    public void testBuildConsensusSpectrum() throws Exception {
//        // iterate over all clusters
//        int index = 0;
//        for (ConsensusSpectraItems cluster : consensusSpectraItems) {
//            ISpectrum consensusSpectrum = cluster.getConcensus();
//            List<ISpectrum> spectra = cluster.getSpectra();
//            if (spectra.size() <= 1)
//                continue;
//
//
//            ISimilarityChecker oldSimilarity = new FrankEtAlDotProductOld();
//            ISimilarityChecker newSimilarity = new FrankEtAlDotProduct();
//
//            for (int index1 = 0; index1 < spectra.size(); index1++) {
//                for (int index2 = index1 + 1; index2 < spectra.size(); index2++) {
//
//                    double oldDotP = oldSimilarity.assessSimilarity(spectra.get(index1), spectra.get(index2));
//                    double newDotP = newSimilarity.assessSimilarity(spectra.get(index1), spectra.get(index2));
//
//                    Assert.assertEquals(oldDotP, newDotP, 0.00001);
//                }
//            }
//
//        }

    }
}






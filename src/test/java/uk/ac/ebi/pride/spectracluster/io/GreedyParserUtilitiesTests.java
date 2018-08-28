package uk.ac.ebi.pride.spectracluster.io;

import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;
import uk.ac.ebi.pride.spectracluster.cluster.GreedySpectralCluster;
import uk.ac.ebi.pride.spectracluster.cluster.ICluster;
import uk.ac.ebi.pride.spectracluster.consensus.BinnedGreedyConsensusSpectrum;
import uk.ac.ebi.pride.spectracluster.similarity.FrankEtAlDotProduct;
import uk.ac.ebi.pride.spectracluster.similarity.ISimilarityChecker;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.util.Defaults;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;

/**
 * Created by jg on 12.05.15.
 */
public class GreedyParserUtilitiesTests {
    private List<ISpectrum> testSpectra;

    @Before
    public void setUp() throws Exception {
        Defaults.resetDefaults();
        File testFile = new File(GreedyParserUtilitiesTests.class.getClassLoader().getResource("spectra_400.0_4.0.mgf").toURI());
        testSpectra = new ArrayList<>();
        ISpectrum[] readSpectra = ParserUtilities.readMGFScans(testFile);

        testSpectra.addAll(Arrays.asList(readSpectra));

        testSpectra.sort(new SpectrumMzComparator());
    }

    @Test
    public void testWriteCluster() throws Exception {
        GreedySpectralCluster cluster = new GreedySpectralCluster("someId");

        for (int i = 0; i < 3; i++)
            cluster.addSpectra(testSpectra.get(i));

        cluster.saveComparisonResult("test_id", 0.7F);

        StringBuilder cgfString = new StringBuilder();

        CGFClusterAppender.INSTANCE.appendCluster(cgfString, cluster);
        String completeString = cgfString.toString();
        String[] lines = completeString.split("\n");

        Assert.assertEquals("BEGIN CLUSTER Id=someId Charge=2 ContainsPeaklist=false", lines[0]);
        Assert.assertEquals("Properties=", lines[1]);
        Assert.assertEquals("ComparisonMatches=0.7:test_id", lines[2]);
        Assert.assertEquals("BEGIN CONSENSUS id=someId class=uk.ac.ebi.pride.spectracluster.consensus.BinnedGreedyConsensusSpectrum nSpec=3 SumCharge=6 SumIntens=0.0 SumMz=1200.8599853515625",
                lines[3]);
    }

    @Test
    public void testReadCluster() throws Exception {
        GreedySpectralCluster cluster = new GreedySpectralCluster("someId");

        for (int i = 0; i < 3; i++)
            cluster.addSpectra(testSpectra.get(i));

        cluster.saveComparisonResult("test_id", 0.7F);

        StringBuilder cgfString = new StringBuilder();

        CGFClusterAppender.INSTANCE.appendCluster(cgfString, cluster);
        ICluster recoveredCluster = ParserUtilities.readSpectralCluster(new LineNumberReader(new StringReader(cgfString.toString())), null);

        Assert.assertTrue(GreedySpectralCluster.class.isInstance(recoveredCluster));
        Assert.assertTrue(BinnedGreedyConsensusSpectrum.class.isInstance(recoveredCluster.getConsensusSpectrumBuilder()));

        Assert.assertEquals("someId", recoveredCluster.getId());
        Assert.assertTrue(recoveredCluster.isKnownComparisonMatch("test_id"));
        Assert.assertFalse(recoveredCluster.isKnownComparisonMatch("test_2"));

        // make sure the consensus spectra are similar
        ISimilarityChecker similarityChecker = new FrankEtAlDotProduct(0.5F);
        double dot = similarityChecker.assessSimilarity(cluster.getConsensusSpectrum(), recoveredCluster.getConsensusSpectrum());
        Assert.assertEquals(1.0, dot, 0.00001);
    }

    @Test
    public void testBinaryReadWriteCluster() throws Exception {
        GreedySpectralCluster cluster = new GreedySpectralCluster("someId");

        for (int i = 0; i < 3; i++)
            cluster.addSpectra(testSpectra.get(i));

        cluster.saveComparisonResult("test_id", 0.7F);

        ByteArrayOutputStream outputStream = new ByteArrayOutputStream();
        ObjectOutputStream objectOutputStream = new ObjectOutputStream(outputStream);

        BinaryClusterAppender.INSTANCE.appendCluster(objectOutputStream, cluster);
        BinaryClusterAppender.INSTANCE.appendEnd(objectOutputStream);
        objectOutputStream.close();
        outputStream.close();

        BinaryClusterIterable binaryClusterIterable = new BinaryClusterIterable(new ObjectInputStream(new ByteArrayInputStream(outputStream.toByteArray())));

        int nClusters = 0;
        ICluster recoveredCluster = null;

        for (ICluster c : binaryClusterIterable) {
            nClusters++;
            recoveredCluster = c;
        }

        Assert.assertEquals(1, nClusters);

        Assert.assertTrue(GreedySpectralCluster.class.isInstance(recoveredCluster));
        Assert.assertTrue(BinnedGreedyConsensusSpectrum.class.isInstance(recoveredCluster.getConsensusSpectrumBuilder()));

        Assert.assertEquals("someId", recoveredCluster.getId());
        Assert.assertTrue(recoveredCluster.isKnownComparisonMatch("test_id"));
        Assert.assertFalse(recoveredCluster.isKnownComparisonMatch("test_2"));

        // make sure the consensus spectra are similar
        ISimilarityChecker similarityChecker = new FrankEtAlDotProduct(0.5F);
        double dot = similarityChecker.assessSimilarity(cluster.getConsensusSpectrum(), recoveredCluster.getConsensusSpectrum());
        Assert.assertEquals(1.0, dot, 0.00001);
    }

    public class SpectrumMzComparator implements Comparator<ISpectrum> {
        @Override
        public int compare(ISpectrum o1, ISpectrum o2) {
            return Float.compare(o1.getPrecursorMz(), o2.getPrecursorMz());
        }
    }

}

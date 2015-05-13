package uk.ac.ebi.pride.spectracluster.io;

import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;
import uk.ac.ebi.pride.spectracluster.cluster.GreedySpectralCluster;
import uk.ac.ebi.pride.spectracluster.cluster.ICluster;
import uk.ac.ebi.pride.spectracluster.consensus.GreedyConsensusSpectrum;
import uk.ac.ebi.pride.spectracluster.engine.GreedySpectralClusteringEngineTest;
import uk.ac.ebi.pride.spectracluster.similarity.FrankEtAlDotProduct;
import uk.ac.ebi.pride.spectracluster.similarity.ISimilarityChecker;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;

import java.io.File;
import java.io.LineNumberReader;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

/**
 * Created by jg on 12.05.15.
 */
public class GreedyParserUtilitiesTests {
    private List<ISpectrum> testSpectra;

    @Before
    public void setUp() throws Exception {
        File testFile = new File(GreedyParserUtilitiesTests.class.getClassLoader().getResource("spectra_400.0_4.0.mgf").toURI());
        testSpectra = new ArrayList<ISpectrum>();
        ISpectrum[] readSpectra = ParserUtilities.readMGFScans(testFile);

        for (ISpectrum s : readSpectra)
            testSpectra.add(s);

        Collections.sort(testSpectra, new SpectrumMzComparator());
    }

    @Test
    public void testWriteCluster() throws Exception {
        GreedySpectralCluster cluster = new GreedySpectralCluster("someId");

        for (int i = 0; i < 3; i++)
            cluster.addSpectra(testSpectra.get(i));

        cluster.saveComparisonResult("test_id", 0.7F);

        StringBuilder cgfString = new StringBuilder();

        CGFClusterAppender.INSTANCE.appendCluster(cgfString, cluster);

        Assert.assertEquals("BEGIN CLUSTER Id=someId Charge=2 ContainsPeaklist=false\n" +
                "ComparisonMatches=0.7:test_id\n" +
                "BEGIN CONSENSUS id=null class=uk.ac.ebi.pride.spectracluster.consensus.GreedyConsensusSpectrum nSpec=3 SumCharge=6 SumIntens=0.0 SumMz=1200.8599853515625\n" +
                "141.871\t60.700\t2\n" +
                "175.087\t502.424\t2\n" +
                "183.206\t90.593\t1\n" +
                "187.017\t83.047\t2\n" +
                "199.081\t6242.819\t3\n" +
                "227.081\t7626.957\t1\n" +
                "251.115\t505.949\t1\n" +
                "271.126\t275.809\t3\n" +
                "274.202\t1892.336\t3\n" +
                "284.153\t567.234\t1\n" +
                "339.931\t1165.851\t4\n" +
                "345.361\t2630.145\t1\n" +
                "381.562\t941.910\t3\n" +
                "390.864\t3245.915\t2\n" +
                "391.964\t730.179\t4\n" +
                "400.778\t149.172\t2\n" +
                "426.472\t386.397\t1\n" +
                "454.204\t2156.580\t1\n" +
                "458.489\t411.352\t1\n" +
                "497.472\t190.703\t1\n" +
                "514.853\t670.884\t1\n" +
                "522.169\t362.603\t1\n" +
                "525.263\t610.760\t1\n" +
                "548.516\t668.679\t1\n" +
                "572.480\t27566.615\t1\n" +
                "624.380\t372.585\t1\n" +
                "661.562\t320.121\t1\n" +
                "669.747\t101.410\t1\n" +
                "684.917\t144.507\t1\n" +
                "685.598\t1301.144\t1\n" +
                "700.317\t2.702\t1\n" +
                "715.665\t23.783\t1\n" +
                "726.909\t4.838\t1\n" +
                "728.840\t70.257\t1\n" +
                "1064.927\t6.158\t1\n" +
                "1150.095\t2.757\t1\n" +
                "1347.983\t3.001\t1\n" +
                "1369.675\t2.447\t1\n" +
                "1392.849\t1.800\t1\n" +
                "1533.479\t4.316\t1\n" +
                "1545.883\t3.463\t1\n" +
                "1607.720\t4.174\t1\n" +
                "END CONSENSUS\n" +
                "BEGIN IONS\n" +
                "TITLE=id=1247848\n" +
                "PEPMASS=400.25\n" +
                "CHARGE=2.0+\n" +
                "USER00=1247848\n" +
                "SEQ=LLGGLAVR\n" +
                "END IONS\n" +
                "BEGIN IONS\n" +
                "TITLE=id=44905\n" +
                "PEPMASS=400.29998779296875\n" +
                "CHARGE=2.0+\n" +
                "USER00=44905\n" +
                "SEQ=KLLQAAR\n" +
                "END IONS\n" +
                "BEGIN IONS\n" +
                "TITLE=id=74575\n" +
                "PEPMASS=400.30999755859375\n" +
                "CHARGE=2.0+\n" +
                "USER00=74575\n" +
                "SEQ=KVIPKSK\n" +
                "END IONS\n" +
                "END CLUSTER\n", cgfString.toString());
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
        Assert.assertTrue(GreedyConsensusSpectrum.class.isInstance(recoveredCluster.getConsensusSpectrumBuilder()));

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

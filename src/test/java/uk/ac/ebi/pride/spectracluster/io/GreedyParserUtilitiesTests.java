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
                "126.038\t1.773\t1\n" +
                "128.970\t18.712\t2\n" +
                "130.121\t7.151\t1\n" +
                "136.021\t36.830\t1\n" +
                "137.982\t39.130\t1\n" +
                "141.871\t38.104\t2\n" +
                "153.995\t33.330\t1\n" +
                "155.031\t28.910\t1\n" +
                "155.795\t2.837\t1\n" +
                "157.094\t10.139\t1\n" +
                "158.138\t11.793\t1\n" +
                "159.172\t41.200\t1\n" +
                "165.149\t14.550\t1\n" +
                "173.089\t4.011\t1\n" +
                "175.087\t315.394\t2\n" +
                "183.206\t78.050\t1\n" +
                "187.017\t52.132\t2\n" +
                "190.058\t9.414\t1\n" +
                "191.151\t30.730\t1\n" +
                "197.296\t20.268\t1\n" +
                "199.081\t2448.164\t3\n" +
                "200.045\t40.620\t1\n" +
                "203.490\t20.230\t1\n" +
                "214.953\t37.330\t1\n" +
                "220.018\t4.408\t1\n" +
                "223.608\t113.512\t2\n" +
                "225.379\t9.064\t1\n" +
                "227.081\t6571.000\t1\n" +
                "228.048\t85.985\t3\n" +
                "229.366\t21.990\t1\n" +
                "232.962\t90.230\t1\n" +
                "234.093\t41.060\t1\n" +
                "234.980\t2.031\t1\n" +
                "240.222\t46.610\t1\n" +
                "241.242\t5.416\t1\n" +
                "243.050\t5.611\t1\n" +
                "249.081\t11.596\t2\n" +
                "251.115\t435.900\t1\n" +
                "251.920\t44.920\t2\n" +
                "255.390\t69.010\t1\n" +
                "259.167\t115.000\t1\n" +
                "269.239\t45.160\t1\n" +
                "271.126\t108.160\t3\n" +
                "273.125\t12.182\t1\n" +
                "274.202\t742.093\t3\n" +
                "279.021\t5.353\t1\n" +
                "279.994\t23.763\t1\n" +
                "281.313\t41.940\t1\n" +
                "282.309\t4.093\t1\n" +
                "283.516\t4.253\t1\n" +
                "284.153\t488.700\t1\n" +
                "286.229\t13.490\t1\n" +
                "287.069\t132.400\t1\n" +
                "288.384\t155.400\t1\n" +
                "290.571\t2.326\t1\n" +
                "292.268\t56.570\t1\n" +
                "293.405\t12.375\t1\n" +
                "296.203\t100.800\t1\n" +
                "297.079\t105.400\t1\n" +
                "299.188\t104.128\t2\n" +
                "300.337\t42.408\t3\n" +
                "301.503\t18.117\t1\n" +
                "304.874\t155.661\t2\n" +
                "306.266\t2.128\t1\n" +
                "306.870\t2.229\t1\n" +
                "309.838\t2.580\t1\n" +
                "311.063\t8.450\t1\n" +
                "313.418\t49.250\t1\n" +
                "314.435\t16.890\t1\n" +
                "315.163\t11.072\t1\n" +
                "317.035\t6.315\t1\n" +
                "318.034\t9.288\t1\n" +
                "318.955\t10.293\t1\n" +
                "322.768\t13.449\t1\n" +
                "324.096\t18.151\t1\n" +
                "324.796\t7.355\t1\n" +
                "325.764\t54.529\t2\n" +
                "327.336\t21.986\t2\n" +
                "328.945\t80.277\t3\n" +
                "331.032\t5.227\t1\n" +
                "336.077\t2.431\t1\n" +
                "339.025\t62.016\t2\n" +
                "339.931\t264.471\t4\n" +
                "341.233\t288.575\t2\n" +
                "342.032\t49.363\t2\n" +
                "342.671\t10.767\t1\n" +
                "344.605\t4.832\t1\n" +
                "345.361\t2266.000\t1\n" +
                "346.545\t14.337\t1\n" +
                "348.374\t15.299\t2\n" +
                "349.111\t7.604\t1\n" +
                "355.209\t161.982\t2\n" +
                "356.064\t13.075\t1\n" +
                "356.740\t115.022\t2\n" +
                "357.729\t40.753\t2\n" +
                "359.001\t14.575\t2\n" +
                "360.271\t20.040\t1\n" +
                "361.722\t5.759\t1\n" +
                "363.524\t22.768\t1\n" +
                "364.545\t41.770\t1\n" +
                "365.070\t88.304\t3\n" +
                "365.756\t2.122\t1\n" +
                "366.256\t63.309\t1\n" +
                "369.402\t79.560\t1\n" +
                "370.505\t21.020\t1\n" +
                "371.223\t11.440\t1\n" +
                "372.072\t33.010\t1\n" +
                "372.807\t4.293\t1\n" +
                "377.865\t167.200\t1\n" +
                "379.823\t73.153\t3\n" +
                "380.599\t29.467\t2\n" +
                "381.562\t369.376\t3\n" +
                "382.375\t225.183\t2\n" +
                "383.259\t91.594\t3\n" +
                "384.155\t23.198\t2\n" +
                "385.132\t250.220\t3\n" +
                "386.321\t154.000\t1\n" +
                "387.763\t84.300\t1\n" +
                "390.864\t2037.606\t2\n" +
                "391.964\t165.640\t4\n" +
                "393.305\t164.626\t2\n" +
                "397.797\t6.539\t1\n" +
                "399.893\t21.356\t1\n" +
                "400.778\t93.642\t2\n" +
                "401.309\t4.194\t1\n" +
                "402.029\t21.429\t1\n" +
                "403.009\t24.731\t2\n" +
                "409.498\t33.200\t1\n" +
                "418.595\t2.485\t1\n" +
                "426.472\t332.900\t1\n" +
                "429.581\t107.000\t1\n" +
                "441.320\t29.990\t1\n" +
                "442.307\t56.550\t1\n" +
                "451.534\t3.049\t1\n" +
                "452.142\t21.700\t1\n" +
                "454.204\t1858.000\t1\n" +
                "458.489\t354.400\t1\n" +
                "459.321\t2.248\t1\n" +
                "462.538\t99.180\t1\n" +
                "466.363\t85.720\t1\n" +
                "474.351\t13.405\t1\n" +
                "476.091\t70.120\t1\n" +
                "483.052\t67.380\t1\n" +
                "484.444\t3.177\t1\n" +
                "485.412\t11.922\t1\n" +
                "486.893\t9.758\t1\n" +
                "488.589\t41.840\t1\n" +
                "490.984\t7.394\t1\n" +
                "497.472\t164.300\t1\n" +
                "507.416\t46.107\t2\n" +
                "510.305\t35.930\t1\n" +
                "511.265\t32.280\t1\n" +
                "512.712\t83.080\t1\n" +
                "514.853\t578.000\t1\n" +
                "515.604\t303.000\t1\n" +
                "519.411\t6.953\t1\n" +
                "522.169\t312.400\t1\n" +
                "525.263\t526.200\t1\n" +
                "527.344\t28.460\t1\n" +
                "532.219\t3.621\t1\n" +
                "532.975\t34.690\t1\n" +
                "537.811\t73.890\t1\n" +
                "538.517\t47.880\t1\n" +
                "541.026\t3.437\t1\n" +
                "542.518\t155.774\t2\n" +
                "543.738\t33.500\t1\n" +
                "546.559\t22.100\t1\n" +
                "548.516\t576.100\t1\n" +
                "551.160\t3.396\t1\n" +
                "554.365\t64.860\t1\n" +
                "558.537\t6.396\t1\n" +
                "570.724\t40.520\t1\n" +
                "572.480\t23750.000\t1\n" +
                "573.494\t225.200\t1\n" +
                "578.451\t41.110\t1\n" +
                "581.586\t3.695\t1\n" +
                "583.102\t39.200\t1\n" +
                "583.653\t4.648\t1\n" +
                "584.993\t4.281\t1\n" +
                "596.283\t3.297\t1\n" +
                "597.379\t77.860\t1\n" +
                "615.377\t2.863\t1\n" +
                "616.012\t56.870\t1\n" +
                "617.550\t4.300\t1\n" +
                "624.380\t321.000\t1\n" +
                "628.459\t10.710\t1\n" +
                "634.858\t1.986\t1\n" +
                "637.407\t52.960\t1\n" +
                "661.562\t275.800\t1\n" +
                "665.370\t43.070\t1\n" +
                "669.747\t87.370\t1\n" +
                "671.591\t6.597\t1\n" +
                "684.917\t124.500\t1\n" +
                "685.598\t1121.000\t1\n" +
                "695.650\t3.290\t1\n" +
                "700.317\t2.328\t1\n" +
                "715.665\t20.490\t1\n" +
                "726.909\t4.169\t1\n" +
                "728.840\t60.530\t1\n" +
                "1064.927\t5.305\t1\n" +
                "1150.095\t2.375\t1\n" +
                "1347.983\t2.586\t1\n" +
                "1369.675\t2.108\t1\n" +
                "1392.849\t1.551\t1\n" +
                "1533.479\t3.718\t1\n" +
                "1545.883\t2.983\t1\n" +
                "1607.720\t3.597\t1\n" +
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

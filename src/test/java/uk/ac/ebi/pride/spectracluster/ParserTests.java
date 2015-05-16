package uk.ac.ebi.pride.spectracluster;

import junit.framework.*;
import org.junit.*;
import org.junit.Assert;
import org.junit.Test;
import uk.ac.ebi.pride.spectracluster.cluster.ICluster;
import uk.ac.ebi.pride.spectracluster.io.CGFClusterAppender;
import uk.ac.ebi.pride.spectracluster.io.CGFSpectrumIterable;
import uk.ac.ebi.pride.spectracluster.spectrum.*;
import uk.ac.ebi.pride.spectracluster.util.ClusteringTestUtilities;
import uk.ac.ebi.pride.spectracluster.io.MGFSpectrumIterable;
import uk.ac.ebi.pride.spectracluster.io.ParserUtilities;

import java.io.File;
import java.io.LineNumberReader;
import java.io.StringReader;
import java.net.URISyntaxException;
import java.net.URL;
import java.util.*;

/**
 * uk.ac.ebi.pride.spectracluster.ParserTests
 *
 * @author Steve Lewis
 * @date 5/12/13
 */
public class ParserTests {
    public static final String MGF_RESOURCE = "res://spectra_400.0_4.0.mgf";

    public static final String[] PEPMASS_STRINGS = {
            "PEPMASS=459.17000000000002 8795.7734375",
            "PEPMASS=445.20999999999998 7451.50732421875",
            "PEPMASS=476.27999999999997 5728.51513671875",
    };
    public static final double[] PEPMASSES = {
            459.17000000000002,
            445.20999999999998,
            476.27999999999997,
    };

    @Test
    public void testPepMassParse() {
        for (int i = 0; i < PEPMASS_STRINGS.length; i++) {
            double mass = ParserUtilities.parsePepMassLine(PEPMASS_STRINGS[i]);
            Assert.assertEquals(PEPMASSES[i], mass, 0.0001);
        }
    }

    @Test
    public void testScansRead() {
        LineNumberReader inp = ParserUtilities.getDescribedReader(MGF_RESOURCE);
        Assert.assertNotNull(inp);  // make sure it exists
        ISpectrum[] spectralClusters = ParserUtilities.readMGFScans(inp);


        Assert.assertEquals(210, spectralClusters.length);

    }


    @Test
    public void testMGFRead() {
        LineNumberReader inp = ParserUtilities.getDescribedReader(MGF_RESOURCE);
        Assert.assertNotNull(inp);  // make sure it exists

        ISpectrum spectralCluster = ParserUtilities.readMGFScan(inp);
        String id;
        int nCluster = 0;
        int[] ids = {1247848, 44905, 74575};

        while (spectralCluster != null) {
            id = spectralCluster.getId();
            Assert.assertNotNull(id);  // make sure it exists
            if (nCluster < ids.length) {
                Assert.assertEquals(ids[nCluster], Integer.parseInt(id));
            }
            String label = spectralCluster.getId();
            Assert.assertNotNull(label);  // make sure it exists

            final List<IPeak> peaks1 = spectralCluster.getPeaks();
            IPeak[] peaks = peaks1.toArray(new IPeak[peaks1.size()]);
            if (peaks.length < 2)
                Assert.assertTrue(peaks.length > 2);  // make there are some peaks
            spectralCluster = ParserUtilities.readMGFScan(inp);
            nCluster++;
        }

        //noinspection UnusedDeclaration,UnusedAssignment
        spectralCluster = ParserUtilities.readMGFScan(inp);
        Assert.assertEquals(210, nCluster);
    }

    /**
     * test reading an MGF Scan and that all ids are unique
     */
    //  @Test
    @SuppressWarnings("UnusedDeclaration")
    public void testUniqueIds() {
        LineNumberReader is = ParserUtilities.getDescribedReader(MGF_RESOURCE);
        testMGFStream(is);

    }

    /**
     * test reading an MGF Scan and that all ids are unique
     */
    protected void testMGFStream(final LineNumberReader inp) {
        Set<String> seenIds = new HashSet<String>();
        ISpectrum ISpectralCluster = ParserUtilities.readMGFScan(inp);
        String id;

        while (ISpectralCluster != null) {

            id = ISpectralCluster.getId();
            Assert.assertTrue(!seenIds.contains(id));   // make sure is is unique
            seenIds.add(id);
            Assert.assertNotNull(id);  // make sure it exists
            String label = ISpectralCluster.getId();
            Assert.assertNotNull(label);  // make sure it exists

            final List<IPeak> peaks1 = ISpectralCluster.getPeaks();
            IPeak[] peaks = peaks1.toArray(new IPeak[peaks1.size()]);
            if (peaks.length < 2)
                Assert.assertTrue(peaks.length > 2);  // make there are some peaks

            ISpectralCluster = ParserUtilities.readMGFScan(inp);
        }

    }


    public static final String CLUSTER_STRING =
            "BEGIN CLUSTER Id=1234 Charge=2 ContainsPeaklist=true\n" +
                    "BEGIN IONS\n" +
                    "TITLE=111749206\n" +
                    "PEPMASS=992.969076\n" +
                    "CHARGE=2+\n" +
                    "320.38278\t3.5\n" +
                    "433.50967\t7.1\n" +
                    "548.22388\t13.5\n" +
                    "685.34204\t46.8\n" +
                    "719.39966\t28.8\n" +
                    "800.37585\t116.0\n" +
                    "914.90295\t598.3\n" +
                    "1038.13965\t75.1\n" +
                    "1092.30835\t57.6\n" +
                    "1193.39893\t63.4\n" +
                    "1340.67883\t10.0\n" +
                    "1437.63416\t79.5\n" +
                    "1552.65173\t51.5\n" +
                    "1640.57629\t18.0\n" +
                    "1738.20154\t2.2\n" +
                    "1810.59561\t73.2\n" +
                    "351.49768\t3.2\n" +
                    "392.05184\t5.4\n" +
                    "555.98425\t5.2\n" +
                    "658.36835\t4.5\n" +
                    "776.53467\t18.2\n" +
                    "833.44678\t21.3\n" +
                    "929.47546\t67.6\n" +
                    "1029.40173\t28.0\n" +
                    "1185.19641\t41.4\n" +
                    "1226.62695\t6.4\n" +
                    "1300.23047\t5.5\n" +
                    "1419.55884\t12.8\n" +
                    "1493.34607\t7.6\n" +
                    "1654.59021\t17.2\n" +
                    "1773.66199\t2.2\n" +
                    "1792.51331\t5.5\n" +
                    "345.3475\t2.7\n" +
                    "446.2587\t3.5\n" +
                    "537.62683\t4.2\n" +
                    "638.17078\t3.3\n" +
                    "768.0047\t6.6\n" +
                    "792.46338\t17.8\n" +
                    "954.98804\t20.3\n" +
                    "1074.37646\t12.7\n" +
                    "1144.33741\t11.9\n" +
                    "1215.49731\t5.2\n" +
                    "1355.49084\t4.9\n" +
                    "1427.69409\t7.7\n" +
                    "1570.2439\t5.1\n" +
                    "1624.04639\t1.2\n" +
                    "1711.16919\t1.7\n" +
                    "1867.42566\t3.5\n" +
                    "326.44226\t2.3\n" +
                    "474.4689\t2.0\n" +
                    "517.2478\t3.9\n" +
                    "677.99765\t3.2\n" +
                    "702.99445\t6.3\n" +
                    "878.41162\t15.1\n" +
                    "893.47498\t12.2\n" +
                    "1022.74512\t6.5\n" +
                    "1124.28613\t7.7\n" +
                    "1233.06348\t2.7\n" +
                    "1368.21655\t3.6\n" +
                    "1443.50793\t5.6\n" +
                    "1534.57971\t3.1\n" +
                    "359.19135\t2.1\n" +
                    "399.28809\t0.8\n" +
                    "563.25531\t3.9\n" +
                    "648.01117\t2.3\n" +
                    "713.17712\t4.3\n" +
                    "807.7666\t6.7\n" +
                    "937.48035\t11.8\n" +
                    "1056.83716\t5.3\n" +
                    "1150.82617\t7.2\n" +
                    "1210.39685\t2.5\n" +
                    "1375.59143\t1.4\n" +
                    "1405.89978\t5.5\n" +
                    "1580.08032\t0.9\n" +
                    "363.18628\t1.8\n" +
                    "439.0647\t0.5\n" +
                    "569.07227\t2.4\n" +
                    "614.20618\t1.9\n" +
                    "784.67834\t4.2\n" +
                    "872.84119\t5.6\n" +
                    "961.63953\t10.6\n" +
                    "1010.22925\t3.4\n" +
                    "1276.5459\t0.8\n" +
                    "1473.61816\t3.8\n" +
                    "377.29718\t1.4\n" +
                    "522.4613\t1.6\n" +
                    "598.22717\t1.5\n" +
                    "737.16754\t3.7\n" +
                    "854.08643\t4.8\n" +
                    "951.91309\t9.9\n" +
                    "1013.76221\t2.3\n" +
                    "1282.32336\t0.4\n" +
                    "1463.65723\t2.7\n" +
                    "286.06775\t1.1\n" +
                    "529.05591\t1.4\n" +
                    "668.30383\t1.3\n" +
                    "695.40674\t2.3\n" +
                    "826.54065\t4.5\n" +
                    "905.2561\t9.3\n" +
                    "1045.7345\t1.6\n" +
                    "1395.59705\t2.3\n" +
                    "632.89044\t1.1\n" +
                    "813.09265\t4.0\n" +
                    "944.09473\t8.7\n" +
                    "1484.94495\t1.7\n" +
                    "848.66418\t3.4\n" +
                    "926.51599\t7.7\n" +
                    "820.80176\t1.5\n" +
                    "891.08875\t2.2\n" +
                    "896.86523\t4.7\n" +
                    "918.5918\t5.3\n" +
                    "921.96448\t4.8\n" +
                    "END IONS\n" +
                    "\n" +
                    "BEGIN IONS\n" +
                    "TITLE=111835642\n" +
                    "PEPMASS=836.770392\n" +
                    "CHARGE=2+\n" +
                    "314.9687\t24.1\n" +
                    "367.22729\t63.7\n" +
                    "525.2395\t136.4\n" +
                    "630.22461\t80.6\n" +
                    "639.62586\t204.4\n" +
                    "779.00055\t691.7\n" +
                    "914.20454\t47.9\n" +
                    "1025.3252\t104.6\n" +
                    "1039.37844\t194.2\n" +
                    "1208.85471\t235.8\n" +
                    "1277.31628\t75.6\n" +
                    "1357.55396\t23.3\n" +
                    "1466.87607\t184.5\n" +
                    "1543.42493\t11.1\n" +
                    "296.23474\t16.4\n" +
                    "383.32849\t55.6\n" +
                    "454.12506\t49.6\n" +
                    "605.03215\t60.8\n" +
                    "648.16364\t121.9\n" +
                    "770.06268\t272.6\n" +
                    "927.25006\t47.3\n" +
                    "1012.36598\t71.0\n" +
                    "1045.23138\t27.2\n" +
                    "1148.41296\t69.7\n" +
                    "1329.29553\t54.5\n" +
                    "1365.14081\t20.9\n" +
                    "1448.31421\t29.4\n" +
                    "1585.5741\t2.8\n" +
                    "266.03751\t12.8\n" +
                    "349.53122\t26.0\n" +
                    "471.02705\t34.1\n" +
                    "634.69462\t58.7\n" +
                    "675.52144\t94.2\n" +
                    "788.31941\t238.0\n" +
                    "1008.14008\t49.5\n" +
                    "1103.40063\t21.9\n" +
                    "1154.08516\t36.1\n" +
                    "1283.24402\t45.6\n" +
                    "1370.29053\t19.5\n" +
                    "1462.49304\t22.1\n" +
                    "236.31824\t12.5\n" +
                    "352.70169\t15.4\n" +
                    "465.15979\t33.9\n" +
                    "606.28113\t45.6\n" +
                    "726.64688\t90.9\n" +
                    "775.38641\t228.8\n" +
                    "976.18622\t48.6\n" +
                    "1112.54065\t21.9\n" +
                    "1168.42981\t34.9\n" +
                    "1308.53944\t25.0\n" +
                    "1422.2926\t19.0\n" +
                    "1524.28601\t18.5\n" +
                    "272.37097\t7.1\n" +
                    "364.62837\t13.3\n" +
                    "448.69327\t25.2\n" +
                    "546.30067\t24.0\n" +
                    "674.05347\t81.8\n" +
                    "784.20868\t83.0\n" +
                    "1029.59094\t46.7\n" +
                    "1085.91934\t21.8\n" +
                    "1201.07764\t25.3\n" +
                    "1331.6908\t21.6\n" +
                    "1361.49133\t17.7\n" +
                    "1437.93079\t18.4\n" +
                    "304.23053\t6.2\n" +
                    "359.67707\t9.7\n" +
                    "534.4245\t18.7\n" +
                    "561.73705\t18.9\n" +
                    "645.16083\t21.3\n" +
                    "797.0824\t71.1\n" +
                    "985.86127\t36.1\n" +
                    "1110.21558\t16.3\n" +
                    "1223.90367\t20.1\n" +
                    "1241.79541\t20.0\n" +
                    "1347.42554\t11.5\n" +
                    "1489.31604\t11.8\n" +
                    "319.43243\t6.2\n" +
                    "372.79501\t9.1\n" +
                    "519.28296\t18.3\n" +
                    "594.05975\t18.9\n" +
                    "800.99542\t58.9\n" +
                    "1035.07642\t30.3\n" +
                    "1120.3739\t14.1\n" +
                    "1204.41431\t19.2\n" +
                    "1260.02792\t19.3\n" +
                    "1416.50732\t11.2\n" +
                    "1508.11396\t11.6\n" +
                    "297.91104\t6.1\n" +
                    "389.7259\t8.4\n" +
                    "508.17035\t14.3\n" +
                    "556.32495\t12.4\n" +
                    "804.55815\t54.4\n" +
                    "1021.3587\t30.0\n" +
                    "1128.38923\t13.2\n" +
                    "1191.34019\t17.1\n" +
                    "1272.422\t16.0\n" +
                    "1412.48401\t10.0\n" +
                    "1505.99609\t8.1\n" +
                    "333.31879\t6.1\n" +
                    "356.10223\t7.7\n" +
                    "502.487\t13.4\n" +
                    "569.93613\t9.7\n" +
                    "806.78714\t53.1\n" +
                    "1014.95892\t29.0\n" +
                    "1115.33826\t12.1\n" +
                    "1177.9735\t15.3\n" +
                    "1289.56226\t14.1\n" +
                    "1423.75806\t9.6\n" +
                    "1445.64722\t7.9\n" +
                    "312.14413\t4.2\n" +
                    "436.19591\t7.6\n" +
                    "532.19373\t13.2\n" +
                    "577.7585\t9.7\n" +
                    "809.91817\t46.5\n" +
                    "1032.07983\t28.1\n" +
                    "1160.39307\t11.2\n" +
                    "1324.18494\t13.0\n" +
                    "1404.40076\t8.8\n" +
                    "1472.29492\t7.5\n" +
                    "327.14075\t2.3\n" +
                    "345.11789\t4.8\n" +
                    "374.91437\t5.1\n" +
                    "396.23907\t2.1\n" +
                    "401.77732\t7.2\n" +
                    "439.56729\t7.9\n" +
                    "457.47968\t6.8\n" +
                    "514.9671\t2.9\n" +
                    "529.93518\t11.2\n" +
                    "537.354\t6.9\n" +
                    "543.49097\t6.3\n" +
                    "553.38629\t9.0\n" +
                    "584.47426\t9.2\n" +
                    "596.03503\t8.9\n" +
                    "600.14368\t4.1\n" +
                    "600.7699\t3.5\n" +
                    "793.18739\t41.1\n" +
                    "981.10187\t16.4\n" +
                    "994.91186\t19.1\n" +
                    "1002.33868\t16.2\n" +
                    "1004.53656\t18.4\n" +
                    "1017.19476\t18.6\n" +
                    "1200.24744\t10.8\n" +
                    "1212.46191\t7.6\n" +
                    "1236.0813\t8.2\n" +
                    "1248.35254\t9.2\n" +
                    "1451.9259\t4.6\n" +
                    "1477.57458\t3.8\n" +
                    "1496.48669\t5.8\n" +
                    "1515.86047\t4.1\n" +
                    "1526.26294\t2.2\n" +
                    "1533.75842\t2.4\n" +
                    "END IONS\n" +
                    "END CLUSTER\n";

    /**
     * test reading an MGF Scan and that all ids are unique
     */
    @Test
    public void testClusterParse() throws Exception {
        LineNumberReader is = new LineNumberReader(new StringReader(CLUSTER_STRING));
        ICluster[] scs = ParserUtilities.readSpectralCluster(is);
        is.close();
        Assert.assertEquals(1, scs.length);
        ICluster sc = scs[0];
        Assert.assertEquals("1234", sc.getId());
        Assert.assertEquals(2, sc.getPrecursorCharge(), 0.001);
        Assert.assertEquals(2, sc.getClusteredSpectraCount());
        List<ISpectrum> clusteredSpectra = sc.getClusteredSpectra();
        Assert.assertEquals(sc.getClusteredSpectraCount(), clusteredSpectra.size());
        for (ISpectrum scc : clusteredSpectra) {
            Assert.assertEquals(2, scc.getPrecursorCharge(), 0.001);
            final List<IPeak> peaks1 = scc.getPeaks();
            IPeak[] peaks = peaks1.toArray(new IPeak[peaks1.size()]);
            Assert.assertTrue(peaks.length > 10);
        }
    }

    /**
     * test reading an MGF Scan and that all ids are unique
     */
    @Test
    public void testClusterWrite() throws Exception {
        LineNumberReader is = new LineNumberReader(new StringReader(CLUSTER_STRING));
        ICluster[] scs = ParserUtilities.readSpectralCluster(is);
        is.close();
        Assert.assertEquals(1, scs.length);
        ICluster sc = scs[0];

        StringBuilder sb = new StringBuilder();
        final CGFClusterAppender clusterAppender = CGFClusterAppender.INSTANCE;
        clusterAppender.appendCluster(sb, sc);
        String sc2Str = sb.toString();
        is = new LineNumberReader(new StringReader(sc2Str));
        scs = ParserUtilities.readSpectralCluster(is);
        is.close();
        Assert.assertEquals(1, scs.length);
        ICluster sc2 = scs[0];


        boolean equivalent = sc.equivalent(sc2);
        Assert.assertTrue(equivalent);

    }

    /**
     * test reading an MGF Scan and that all ids are unique
     */
    @Test
    public void testClusterIterator() throws Exception {
        LineNumberReader is = new LineNumberReader(new StringReader(CLUSTER_STRING));
        ICluster[] scs = ParserUtilities.readSpectralCluster(is);
        is.close();
        Assert.assertEquals(1, scs.length);


        is = new LineNumberReader(new StringReader(CLUSTER_STRING));
        CGFSpectrumIterable mgi = new CGFSpectrumIterable(is);
        List<ICluster> holder = new ArrayList<ICluster>();
        for (ICluster sc : mgi) {
            holder.add(sc);
        }
        Assert.assertEquals(1, holder.size());

        ICluster sc = scs[0];
        ICluster sc2 = holder.get(0);

        boolean equivalent = sc.equivalent(sc2);
        Assert.assertTrue(equivalent);

    }

    /**
     * test reading an MGF Scan conventionally and with an iterator
     */
    @Test
    public void testSpectrumIterator() throws Exception {

        List<? extends ISpectrum> spectra = ClusteringTestUtilities.readISpectraFromResource();

        // load a file contains a list of clusters
        URL url;
        url = ParserTests.class.getClassLoader().getResource(ClusteringTestUtilities.SAMPLE_MGF_FILE);
        if (url == null) {
            throw new IllegalStateException("no file for input found!");
        }
        File inputFile;
        try {
            inputFile = new File(url.toURI());
        } catch (URISyntaxException e) {
            throw new RuntimeException(e);

        }

        List<ISpectrum> holder = new ArrayList<ISpectrum>();
        MGFSpectrumIterable itr = new MGFSpectrumIterable(inputFile);
        for (ISpectrum sp : itr) {
            holder.add(sp);
        }
        Assert.assertEquals(holder.size(), spectra.size());

        Collections.sort(holder);
        Collections.sort(spectra);

        // bettrer be the same
        for (int i = 0; i < holder.size(); i++) {
            Assert.assertTrue(holder.get(i).equivalent(spectra.get(i)));

        }

    }

}

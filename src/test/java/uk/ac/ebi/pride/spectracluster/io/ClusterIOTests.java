package uk.ac.ebi.pride.spectracluster.io;

import org.junit.Assert;
import org.junit.Test;
import uk.ac.ebi.pride.spectracluster.cluster.ICluster;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.util.ClusterUtilities;

import java.io.LineNumberReader;
import java.io.StringReader;

/**
 * uk.ac.ebi.pride.spectracluster.io.ClusterIOTests
 *
 * @author Steve Lewis
 * @date 04/06/2014
 */
public class ClusterIOTests {

    @Test
    public void testSpectrumAppender() {
        LineNumberReader rdr = new LineNumberReader(new StringReader(TEST_SPECTRUM));
        final ISpectrum spc = ParserUtilities.readMGFScan(rdr);
        StringBuilder bldr = new StringBuilder();
        MGFSpectrumAppender.INSTANCE.appendSpectrum(bldr, spc);
        final String actual = bldr.toString();
        rdr = new LineNumberReader(new StringReader(actual));
        final ISpectrum spc2 = ParserUtilities.readMGFScan(rdr);
        Assert.assertTrue(spc.equivalent(spc2));
    }


    @Test
    public void testClusterAppender() {
        LineNumberReader rdr = new LineNumberReader(new StringReader(TEST_CLUSTER));
        final ICluster spc = ParserUtilities.readSpectralCluster(rdr)[0];
        StringBuilder bldr = new StringBuilder();
        CGFClusterAppender.INSTANCE.appendCluster(bldr, spc);
        final String actual = bldr.toString();
        rdr = new LineNumberReader(new StringReader(actual));
        final ICluster spc2 = ParserUtilities.readSpectralCluster(rdr)[0];
        Assert.assertTrue(spc.equivalent(spc2));
    }

    @Test
    public void testClusterSpectra() {
        LineNumberReader rdr = new LineNumberReader(new StringReader(TEST_CLUSTER));
        final ICluster cluster = ParserUtilities.readSpectralCluster(rdr)[0];
        String s = ClusterUtilities.mostCommonPeptides(cluster);
        Assert.assertEquals("MLQQGSSSR:1,MSANISETTAM:1,SAKAQNPMR:1,STKAQNPMR:1", s);
    }

    public static final String TEST_SPECTRUM =
            "BEGIN IONS\n" +
                    "TITLE=id=PXD000316;PRIDE_Exp_Complete_Ac_29985.xml;spectrum=139\n" +
                    "PEPMASS=569.2871704101562\n" +
                    "CHARGE=2.0+\n" +
                    "SEQ=STKAQNPMR\n" +
                    "TAXONOMY=4932\n" +
                    "USER02=SwissProt:Q3E757\n" +
                    "175.218\t6.757\n" +
                    "209.261\t1.065\n" +
                    "212.055\t1.615\n" +
                    "213.004\t2.493\n" +
                    "221.216\t1.355\n" +
                    "230.216\t1.406\n" +
                    "240.270\t1.002\n" +
                    "259.057\t2.307\n" +
                    "265.186\t2.267\n" +
                    "273.023\t1.230\n" +
                    "277.169\t1.675\n" +
                    "314.311\t3.357\n" +
                    "318.158\t1.399\n" +
                    "322.104\t1.861\n" +
                    "347.379\t1.354\n" +
                    "348.402\t2.782\n" +
                    "350.442\t3.065\n" +
                    "352.175\t1.171\n" +
                    "355.279\t3.069\n" +
                    "359.093\t1.626\n" +
                    "363.239\t3.573\n" +
                    "369.560\t4.243\n" +
                    "370.327\t9.850\n" +
                    "373.242\t1.442\n" +
                    "375.210\t2.488\n" +
                    "376.785\t3.175\n" +
                    "385.374\t3.868\n" +
                    "388.071\t46.316\n" +
                    "392.991\t1.644\n" +
                    "396.796\t34.374\n" +
                    "401.058\t4.642\n" +
                    "406.053\t18.548\n" +
                    "407.362\t6.339\n" +
                    "410.726\t1.852\n" +
                    "419.194\t182.605\n" +
                    "421.118\t5.898\n" +
                    "435.968\t2.797\n" +
                    "437.027\t7.512\n" +
                    "439.326\t3.019\n" +
                    "441.024\t9.153\n" +
                    "447.077\t2.846\n" +
                    "449.054\t1.482\n" +
                    "454.002\t29.399\n" +
                    "455.049\t8.772\n" +
                    "459.143\t11.783\n" +
                    "460.297\t4.998\n" +
                    "464.161\t3.386\n" +
                    "466.989\t2.699\n" +
                    "469.649\t4.568\n" +
                    "471.148\t6.434\n" +
                    "472.963\t6.475\n" +
                    "476.418\t23.710\n" +
                    "477.235\t81.378\n" +
                    "482.295\t2.116\n" +
                    "484.212\t2.770\n" +
                    "487.459\t3.418\n" +
                    "488.231\t7.134\n" +
                    "489.013\t3.451\n" +
                    "495.142\t4.098\n" +
                    "496.449\t16.850\n" +
                    "504.921\t115.574\n" +
                    "512.119\t6.267\n" +
                    "515.155\t110.034\n" +
                    "516.034\t28.157\n" +
                    "519.738\t9.775\n" +
                    "524.737\t4.341\n" +
                    "526.271\t7.775\n" +
                    "528.046\t5.110\n" +
                    "533.188\t167.922\n" +
                    "536.716\t17.425\n" +
                    "537.513\t31.791\n" +
                    "543.047\t24.742\n" +
                    "546.455\t19.497\n" +
                    "551.196\t335.725\n" +
                    "552.009\t167.355\n" +
                    "560.315\t283.051\n" +
                    "568.821\t5.603\n" +
                    "587.225\t8.654\n" +
                    "590.257\t9.117\n" +
                    "594.476\t5.141\n" +
                    "605.119\t32.814\n" +
                    "617.414\t2.193\n" +
                    "620.896\t2.289\n" +
                    "622.414\t1.488\n" +
                    "635.471\t1.722\n" +
                    "646.403\t3.318\n" +
                    "650.099\t1.256\n" +
                    "661.258\t74.083\n" +
                    "676.324\t2.038\n" +
                    "683.342\t6.440\n" +
                    "684.444\t3.766\n" +
                    "685.533\t2.272\n" +
                    "691.421\t2.429\n" +
                    "702.182\t35.849\n" +
                    "719.336\t130.646\n" +
                    "732.290\t74.468\n" +
                    "757.270\t1.660\n" +
                    "764.629\t3.557\n" +
                    "799.503\t2.411\n" +
                    "816.356\t3.831\n" +
                    "834.567\t4.901\n" +
                    "838.377\t1.866\n" +
                    "843.507\t2.739\n" +
                    "855.131\t1.998\n" +
                    "880.513\t1.750\n" +
                    "907.346\t4.848\n" +
                    "914.262\t3.124\n" +
                    "946.245\t7.659\n" +
                    "963.383\t26.006\n" +
                    "969.070\t1.719\n" +
                    "991.630\t8.864\n" +
                    "1098.353\t9.571\n" +
                    "END IONS\n" +
                    "BEGIN IONS\n";


    public static final String TEST_CLUSTER =
            "BEGIN CLUSTER Id=22706 Charge=2 ContainsPeaklist=true\n" +
                    "BEGIN IONS\n" +
                    "TITLE=id=PXD000316;PRIDE_Exp_Complete_Ac_29985.xml;spectrum=139\n" +
                    "PEPMASS=569.2871704101562\n" +
                    "CHARGE=2.0+\n" +
                    "SEQ=STKAQNPMR\n" +
                    "TAXONOMY=4932\n" +
                    "USER02=SwissProt:Q3E757\n" +
                    "175.218\t6.757\n" +
                    "209.261\t1.065\n" +
                    "212.055\t1.615\n" +
                    "213.004\t2.493\n" +
                    "221.216\t1.355\n" +
                    "230.216\t1.406\n" +
                    "240.270\t1.002\n" +
                    "259.057\t2.307\n" +
                    "265.186\t2.267\n" +
                    "273.023\t1.230\n" +
                    "277.169\t1.675\n" +
                    "314.311\t3.357\n" +
                    "318.158\t1.399\n" +
                    "322.104\t1.861\n" +
                    "347.379\t1.354\n" +
                    "348.402\t2.782\n" +
                    "350.442\t3.065\n" +
                    "352.175\t1.171\n" +
                    "355.279\t3.069\n" +
                    "359.093\t1.626\n" +
                    "363.239\t3.573\n" +
                    "369.560\t4.243\n" +
                    "370.327\t9.850\n" +
                    "373.242\t1.442\n" +
                    "375.210\t2.488\n" +
                    "376.785\t3.175\n" +
                    "385.374\t3.868\n" +
                    "388.071\t46.316\n" +
                    "392.991\t1.644\n" +
                    "396.796\t34.374\n" +
                    "401.058\t4.642\n" +
                    "406.053\t18.548\n" +
                    "407.362\t6.339\n" +
                    "410.726\t1.852\n" +
                    "419.194\t182.605\n" +
                    "421.118\t5.898\n" +
                    "435.968\t2.797\n" +
                    "437.027\t7.512\n" +
                    "439.326\t3.019\n" +
                    "441.024\t9.153\n" +
                    "447.077\t2.846\n" +
                    "449.054\t1.482\n" +
                    "454.002\t29.399\n" +
                    "455.049\t8.772\n" +
                    "459.143\t11.783\n" +
                    "460.297\t4.998\n" +
                    "464.161\t3.386\n" +
                    "466.989\t2.699\n" +
                    "469.649\t4.568\n" +
                    "471.148\t6.434\n" +
                    "472.963\t6.475\n" +
                    "476.418\t23.710\n" +
                    "477.235\t81.378\n" +
                    "482.295\t2.116\n" +
                    "484.212\t2.770\n" +
                    "487.459\t3.418\n" +
                    "488.231\t7.134\n" +
                    "489.013\t3.451\n" +
                    "495.142\t4.098\n" +
                    "496.449\t16.850\n" +
                    "504.921\t115.574\n" +
                    "512.119\t6.267\n" +
                    "515.155\t110.034\n" +
                    "516.034\t28.157\n" +
                    "519.738\t9.775\n" +
                    "524.737\t4.341\n" +
                    "526.271\t7.775\n" +
                    "528.046\t5.110\n" +
                    "533.188\t167.922\n" +
                    "536.716\t17.425\n" +
                    "537.513\t31.791\n" +
                    "543.047\t24.742\n" +
                    "546.455\t19.497\n" +
                    "551.196\t335.725\n" +
                    "552.009\t167.355\n" +
                    "560.315\t283.051\n" +
                    "568.821\t5.603\n" +
                    "587.225\t8.654\n" +
                    "590.257\t9.117\n" +
                    "594.476\t5.141\n" +
                    "605.119\t32.814\n" +
                    "617.414\t2.193\n" +
                    "620.896\t2.289\n" +
                    "622.414\t1.488\n" +
                    "635.471\t1.722\n" +
                    "646.403\t3.318\n" +
                    "650.099\t1.256\n" +
                    "661.258\t74.083\n" +
                    "676.324\t2.038\n" +
                    "683.342\t6.440\n" +
                    "684.444\t3.766\n" +
                    "685.533\t2.272\n" +
                    "691.421\t2.429\n" +
                    "702.182\t35.849\n" +
                    "719.336\t130.646\n" +
                    "732.290\t74.468\n" +
                    "757.270\t1.660\n" +
                    "764.629\t3.557\n" +
                    "799.503\t2.411\n" +
                    "816.356\t3.831\n" +
                    "834.567\t4.901\n" +
                    "838.377\t1.866\n" +
                    "843.507\t2.739\n" +
                    "855.131\t1.998\n" +
                    "880.513\t1.750\n" +
                    "907.346\t4.848\n" +
                    "914.262\t3.124\n" +
                    "946.245\t7.659\n" +
                    "963.383\t26.006\n" +
                    "969.070\t1.719\n" +
                    "991.630\t8.864\n" +
                    "1098.353\t9.571\n" +
                    "END IONS\n" +
                    "BEGIN IONS\n" +
                    "TITLE=id=PXD000316;PRIDE_Exp_Complete_Ac_29985.xml;spectrum=171\n" +
                    "PEPMASS=554.2818603515625\n" +
                    "CHARGE=2.0+\n" +
                    "SEQ=SAKAQNPMR\n" +
                    "TAXONOMY=4932\n" +
                    "USER02=SwissProt:P0C0W9\n" +
                    "175.147\t21.350\n" +
                    "183.125\t8.343\n" +
                    "199.697\t3.147\n" +
                    "200.568\t3.724\n" +
                    "217.948\t11.484\n" +
                    "247.214\t6.088\n" +
                    "265.210\t6.607\n" +
                    "310.897\t4.273\n" +
                    "314.297\t3.108\n" +
                    "322.255\t7.020\n" +
                    "340.167\t4.840\n" +
                    "345.933\t2.488\n" +
                    "358.119\t47.006\n" +
                    "370.168\t5.069\n" +
                    "373.378\t3.570\n" +
                    "375.449\t21.039\n" +
                    "376.242\t90.419\n" +
                    "381.677\t3.726\n" +
                    "387.239\t12.372\n" +
                    "388.355\t5.505\n" +
                    "393.129\t5.127\n" +
                    "401.210\t12.590\n" +
                    "402.364\t3.171\n" +
                    "405.469\t4.056\n" +
                    "410.687\t3.662\n" +
                    "413.762\t7.840\n" +
                    "419.217\t432.852\n" +
                    "426.151\t3.173\n" +
                    "428.646\t13.413\n" +
                    "429.339\t28.501\n" +
                    "433.600\t3.250\n" +
                    "438.552\t9.401\n" +
                    "440.392\t7.096\n" +
                    "447.201\t205.972\n" +
                    "454.384\t82.888\n" +
                    "460.841\t8.290\n" +
                    "461.863\t3.138\n" +
                    "464.256\t3.409\n" +
                    "468.877\t5.325\n" +
                    "470.247\t3.348\n" +
                    "473.180\t3.455\n" +
                    "481.376\t19.815\n" +
                    "489.708\t139.018\n" +
                    "490.479\t67.893\n" +
                    "494.161\t12.601\n" +
                    "495.858\t8.954\n" +
                    "501.397\t14.795\n" +
                    "503.637\t15.812\n" +
                    "505.392\t7.513\n" +
                    "508.849\t5.881\n" +
                    "521.587\t19.647\n" +
                    "522.432\t54.259\n" +
                    "533.266\t388.226\n" +
                    "535.409\t77.727\n" +
                    "536.500\t282.404\n" +
                    "539.665\t20.520\n" +
                    "545.636\t818.571\n" +
                    "557.153\t10.193\n" +
                    "560.240\t10.948\n" +
                    "574.509\t9.900\n" +
                    "575.130\t44.609\n" +
                    "608.045\t2.713\n" +
                    "616.320\t3.100\n" +
                    "631.032\t3.201\n" +
                    "643.473\t10.618\n" +
                    "644.279\t5.936\n" +
                    "661.273\t185.431\n" +
                    "668.489\t6.764\n" +
                    "671.261\t16.488\n" +
                    "672.268\t63.220\n" +
                    "688.802\t46.352\n" +
                    "689.416\t230.713\n" +
                    "700.555\t9.229\n" +
                    "715.227\t26.099\n" +
                    "719.349\t6.552\n" +
                    "732.320\t255.615\n" +
                    "775.451\t3.225\n" +
                    "786.548\t5.761\n" +
                    "858.817\t3.128\n" +
                    "907.532\t7.638\n" +
                    "916.193\t12.199\n" +
                    "933.307\t51.850\n" +
                    "1069.639\t3.524\n" +
                    "END IONS\n" +
                    "BEGIN IONS\n" +
                    "TITLE=id=PXD000316;PRIDE_Exp_Complete_Ac_29985.xml;spectrum=195\n" +
                    "PEPMASS=526.2453002929688\n" +
                    "CHARGE=2.0+\n" +
                    "SEQ=MLQQGSSSR\n" +
                    "TAXONOMY=4932\n" +
                    "USER02=SwissProt:P40095\n" +
                    "162.128\t6.814\n" +
                    "175.088\t19.423\n" +
                    "190.077\t254.050\n" +
                    "242.100\t47.319\n" +
                    "244.197\t8.203\n" +
                    "257.181\t5.935\n" +
                    "262.243\t136.816\n" +
                    "274.288\t5.205\n" +
                    "275.117\t7.192\n" +
                    "303.109\t67.721\n" +
                    "313.165\t13.870\n" +
                    "320.188\t9.730\n" +
                    "331.212\t7.279\n" +
                    "349.219\t322.329\n" +
                    "353.316\t7.552\n" +
                    "357.219\t8.246\n" +
                    "365.127\t19.396\n" +
                    "366.085\t6.690\n" +
                    "366.833\t21.620\n" +
                    "370.129\t42.255\n" +
                    "371.342\t10.050\n" +
                    "375.306\t59.374\n" +
                    "399.938\t17.098\n" +
                    "401.528\t10.229\n" +
                    "402.269\t5.165\n" +
                    "411.899\t10.320\n" +
                    "413.125\t14.746\n" +
                    "413.892\t19.342\n" +
                    "423.067\t82.307\n" +
                    "431.625\t833.286\n" +
                    "436.218\t128.859\n" +
                    "442.292\t5.699\n" +
                    "447.772\t7.325\n" +
                    "457.238\t19.998\n" +
                    "461.826\t8.274\n" +
                    "466.279\t10.948\n" +
                    "472.665\t11.073\n" +
                    "475.298\t19.528\n" +
                    "478.901\t12.738\n" +
                    "485.100\t25.199\n" +
                    "487.332\t11.639\n" +
                    "488.752\t22.895\n" +
                    "493.294\t527.999\n" +
                    "496.845\t16.846\n" +
                    "499.351\t34.098\n" +
                    "508.368\t916.133\n" +
                    "517.387\t972.687\n" +
                    "542.187\t7.317\n" +
                    "559.133\t254.092\n" +
                    "585.471\t11.298\n" +
                    "586.613\t5.618\n" +
                    "603.298\t43.201\n" +
                    "604.249\t15.523\n" +
                    "616.268\t158.024\n" +
                    "621.278\t639.677\n" +
                    "675.060\t20.880\n" +
                    "685.225\t66.761\n" +
                    "703.233\t306.453\n" +
                    "731.331\t16.048\n" +
                    "732.337\t30.974\n" +
                    "749.349\t276.516\n" +
                    "754.174\t25.957\n" +
                    "772.128\t74.646\n" +
                    "789.978\t47.536\n" +
                    "826.540\t24.635\n" +
                    "827.756\t10.873\n" +
                    "845.391\t80.630\n" +
                    "862.442\t883.537\n" +
                    "877.223\t17.377\n" +
                    "END IONS\n" +
                    "BEGIN IONS\n" +
                    "TITLE=id=PXD000316;PRIDE_Exp_Complete_Ac_29985.xml;spectrum=317\n" +
                    "PEPMASS=594.2498168945312\n" +
                    "CHARGE=2.0+\n" +
                    "SEQ=MSANISETTAM\n" +
                    "TAXONOMY=4932\n" +
                    "USER02=SwissProt:P06169\n" +
                    "175.076\t2.658\n" +
                    "206.189\t1.719\n" +
                    "219.039\t4.939\n" +
                    "221.118\t4.543\n" +
                    "237.111\t9.027\n" +
                    "246.270\t10.090\n" +
                    "247.160\t1.567\n" +
                    "258.143\t13.794\n" +
                    "259.325\t1.544\n" +
                    "267.244\t4.796\n" +
                    "281.057\t11.069\n" +
                    "286.096\t1.554\n" +
                    "288.243\t2.472\n" +
                    "291.089\t17.048\n" +
                    "293.273\t1.544\n" +
                    "299.228\t2.116\n" +
                    "300.276\t3.812\n" +
                    "314.277\t3.470\n" +
                    "322.165\t11.197\n" +
                    "326.387\t2.477\n" +
                    "338.019\t10.088\n" +
                    "340.205\t3.037\n" +
                    "343.164\t1.803\n" +
                    "344.066\t2.669\n" +
                    "347.209\t2.172\n" +
                    "355.069\t209.970\n" +
                    "355.977\t10.004\n" +
                    "370.984\t30.562\n" +
                    "372.164\t5.646\n" +
                    "376.147\t5.174\n" +
                    "382.975\t1.439\n" +
                    "392.246\t1.836\n" +
                    "403.327\t2.430\n" +
                    "404.894\t1.478\n" +
                    "417.068\t1.837\n" +
                    "418.368\t2.271\n" +
                    "420.214\t5.822\n" +
                    "421.190\t2.111\n" +
                    "424.340\t3.702\n" +
                    "428.303\t3.245\n" +
                    "430.678\t2.234\n" +
                    "435.173\t2.689\n" +
                    "439.106\t94.379\n" +
                    "448.260\t3.317\n" +
                    "451.813\t4.722\n" +
                    "458.158\t3.048\n" +
                    "459.426\t3.034\n" +
                    "463.042\t16.615\n" +
                    "464.238\t3.843\n" +
                    "468.053\t4.547\n" +
                    "470.100\t2.869\n" +
                    "473.006\t3.283\n" +
                    "478.335\t9.358\n" +
                    "479.260\t4.334\n" +
                    "480.105\t10.125\n" +
                    "484.848\t5.205\n" +
                    "487.997\t2.664\n" +
                    "490.261\t6.704\n" +
                    "494.844\t4.620\n" +
                    "504.210\t6.352\n" +
                    "510.398\t4.787\n" +
                    "511.638\t7.756\n" +
                    "513.273\t4.263\n" +
                    "516.139\t4.132\n" +
                    "517.132\t8.969\n" +
                    "518.407\t4.490\n" +
                    "522.554\t8.602\n" +
                    "528.288\t237.567\n" +
                    "530.436\t25.977\n" +
                    "533.136\t199.491\n" +
                    "537.089\t6.420\n" +
                    "539.333\t15.281\n" +
                    "541.259\t7.066\n" +
                    "544.604\t17.770\n" +
                    "550.875\t8.290\n" +
                    "553.024\t10.648\n" +
                    "554.578\t11.541\n" +
                    "556.800\t13.066\n" +
                    "557.892\t22.923\n" +
                    "562.507\t24.082\n" +
                    "568.197\t27.467\n" +
                    "576.246\t190.732\n" +
                    "577.050\t55.162\n" +
                    "583.616\t54.201\n" +
                    "585.413\t110.438\n" +
                    "594.619\t4.371\n" +
                    "598.646\t2.028\n" +
                    "602.262\t3.473\n" +
                    "620.202\t26.839\n" +
                    "631.452\t4.761\n" +
                    "637.293\t16.218\n" +
                    "643.530\t1.744\n" +
                    "655.143\t189.472\n" +
                    "660.706\t2.385\n" +
                    "684.039\t1.903\n" +
                    "690.923\t2.768\n" +
                    "700.503\t3.407\n" +
                    "731.303\t10.325\n" +
                    "732.340\t20.262\n" +
                    "749.244\t117.076\n" +
                    "768.526\t6.803\n" +
                    "778.644\t2.051\n" +
                    "799.980\t1.754\n" +
                    "817.314\t1.867\n" +
                    "832.351\t6.246\n" +
                    "834.507\t2.991\n" +
                    "850.293\t16.743\n" +
                    "875.012\t1.828\n" +
                    "900.295\t4.394\n" +
                    "923.261\t2.509\n" +
                    "933.505\t14.087\n" +
                    "951.320\t30.790\n" +
                    "1004.463\t24.844\n" +
                    "1022.292\t74.318\n" +
                    "END IONS\n" +
                    "END CLUSTER\n";
}

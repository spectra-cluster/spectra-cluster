package uk.ac.ebi.pride.spectracluster.io;

import junit.framework.Assert;
import org.junit.Before;
import org.junit.Test;
import uk.ac.ebi.pride.spectracluster.cluster.ICluster;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.spectrum.KnownProperties;

import java.io.File;
import java.io.StringReader;
import java.util.List;

/**
 * Created by jg on 01.06.17.
 */
public class PrideMgfParameterTest {
    protected File testMgf;

    @Before
    public void setUp() throws Exception {
        testMgf = new File(
                PrideMgfParameterTest.class.getClassLoader().getResource("pride_export_2017.mgf").toURI());
    }

    @Test
    public void testParseMgf() throws Exception {
        ISpectrum[] spectra = ParserUtilities.readMGFScans(testMgf);

        Assert.assertEquals(312, spectra.length);

        String[] newProperties = {KnownProperties.PSM_FDR_SCORES, KnownProperties.PSM_DECOY_STATUS};

        // make sure all spectra have the new properties
        for (ISpectrum s : spectra) {
            for (String property : newProperties) {
                String value = s.getProperty(property);
                Assert.assertNotNull(value);
            }
        }

        // make sure the first and last one are correct
        Assert.assertEquals("decoy:true;decoy:false;decoy:true",
                spectra[0].getProperty(KnownProperties.PSM_DECOY_STATUS));
        Assert.assertEquals("psm_combined_fdr_score=0.001066603460535672;psm_combined_fdr_score=5.806004411790403E-5;psm_combined_fdr_score=0.001066603460535672",
                spectra[0].getProperty(KnownProperties.PSM_FDR_SCORES));

        Assert.assertEquals("decoy:false",
                spectra[311].getProperty(KnownProperties.PSM_DECOY_STATUS));
        Assert.assertEquals("psm_combined_fdr_score=6.69193471563635E-5",
                spectra[311].getProperty(KnownProperties.PSM_FDR_SCORES));
    }

    @Test
    public void testCgfAppender() throws Exception {
        // get the spectra as clusters
        List<ICluster> clusters = ParserUtilities.readMGFClusters(testMgf);

        Assert.assertEquals(312, clusters.size());

        // convert to CGF
        StringBuilder stringBuilder = new StringBuilder();
        CGFClusterAppender.INSTANCE.appendCluster(stringBuilder, clusters.get(0));

        String testString = stringBuilder.toString();

        String[] lines = testString.split("\n");
        Assert.assertEquals(163, lines.length);
        Assert.assertEquals("USER05=decoy:true;decoy:false;decoy:true", lines[9]);
        Assert.assertEquals("USER06=psm_combined_fdr_score=0.001066603460535672;psm_combined_fdr_score=5.806004411790403E-5;psm_combined_fdr_score=0.001066603460535672",
                lines[10]);

        // read the cluster in again
        CGFSpectrumIterable clusterIterator = new CGFSpectrumIterable(new StringReader(testString));
        Assert.assertTrue(clusterIterator.iterator().hasNext());
        ICluster cgfCluster = clusterIterator.iterator().next();

        Assert.assertEquals(1, cgfCluster.getClusteredSpectraCount());
        Assert.assertEquals("decoy:true;decoy:false;decoy:true",
                cgfCluster.getClusteredSpectra().get(0).getProperty(KnownProperties.PSM_DECOY_STATUS));
        Assert.assertEquals("psm_combined_fdr_score=0.001066603460535672;psm_combined_fdr_score=5.806004411790403E-5;psm_combined_fdr_score=0.001066603460535672",
                cgfCluster.getClusteredSpectra().get(0).getProperty(KnownProperties.PSM_FDR_SCORES));
    }

    @Test
    public void testClusteringAppender() throws Exception {
        // get the spectra as clusters
        List<ICluster> clusters = ParserUtilities.readMGFClusters(testMgf);

        StringBuilder stringBuilder = new StringBuilder();

        DotClusterClusterAppender.INSTANCE.appendCluster(stringBuilder, clusters.get(0));

        String clusteringResult = stringBuilder.toString();

        String[] lines = clusteringResult.split("\n");

        Assert.assertEquals(9, lines.length);

        String[] specFields = lines[8].split("\t");

        Assert.assertEquals(10, specFields.length);
        Assert.assertEquals("{\"DECOY\": \"decoy:true;decoy:false;decoy:true\", \"FDR\": \"psm_combined_fdr_score=0.001066603460535672;psm_combined_fdr_score=5.806004411790403E-5;psm_combined_fdr_score=0.001066603460535672\"}",
                specFields[9]);
    }
}

package uk.ac.ebi.pride.spectracluster.io;

import junit.framework.Assert;
import org.junit.Before;
import org.junit.Test;
import uk.ac.ebi.pride.spectracluster.cluster.ICluster;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.spectrum.KnownProperties;

import javax.xml.bind.annotation.adapters.XmlJavaTypeAdapter;
import java.io.File;
import java.io.StringReader;
import java.util.List;

/**
 * Created by jg on 01.06.17.
 */
public class RetentionTimeTest {
    protected File testMgf;

    @Before
    public void setUp() throws Exception {
        testMgf = new File(
                RetentionTimeTest.class.getClassLoader().getResource("kuester_test.mgf").toURI());
    }

    @Test
    public void testParseMgf() throws Exception {
        ISpectrum[] spectra = ParserUtilities.readMGFScans(testMgf);

        Assert.assertEquals(171, spectra.length);

        // make sure all spectra have retention times
        for (ISpectrum s : spectra) {
            String rt = s.getProperty(KnownProperties.RETENTION_TIME);
            Assert.assertNotNull(rt);
        }

        // make sure the first and last one are correct
        Assert.assertEquals("2293.276122", spectra[0].getProperty(KnownProperties.RETENTION_TIME));
        Assert.assertEquals("3600.959897", spectra[170].getProperty(KnownProperties.RETENTION_TIME));
    }

    @Test
    public void testCgfAppender() throws Exception {
        // get the spectra as clusters
        List<ICluster> clusters = ParserUtilities.readMGFClusters(testMgf);

        Assert.assertEquals(171, clusters.size());

        // make sure the retention time is there
        Assert.assertEquals("2293.276122",
                clusters.get(0).getClusteredSpectra().get(0).getProperty(KnownProperties.RETENTION_TIME));

        // convert to CGF
        ICluster cluster = clusters.get(0);
        Assert.assertEquals("2293.276122", cluster.getClusteredSpectra().get(0).getProperty(KnownProperties.RETENTION_TIME));
        cluster.getClusteredSpectra().get(0).setProperty(KnownProperties.MIN_COMPARISONS, "1234");
        cluster.getClusteredSpectra().get(0).setProperty(KnownProperties.ADDING_SCORE, "17.1");
        StringBuilder stringBuilder = new StringBuilder();
        CGFClusterAppender.INSTANCE.appendCluster(stringBuilder, cluster);

        String testString = stringBuilder.toString();

        String[] lines = testString.split("\n");

        if (lines.length != 163) {
            System.out.print(testString);
        }

        Assert.assertEquals(163, lines.length);
        Assert.assertEquals("RTINSECONDS=2293.276122", lines[8]);
        Assert.assertEquals("MIN_COMP=1234", lines[9]);
        Assert.assertEquals("ADDING_SCORE=17.1", lines[10]);

        // read the cluster in again
        CGFSpectrumIterable clusterIterator = new CGFSpectrumIterable(new StringReader(testString));
        Assert.assertTrue(clusterIterator.iterator().hasNext());
        ICluster cgfCluster = clusterIterator.iterator().next();

        Assert.assertEquals(1, cgfCluster.getClusteredSpectraCount());
        Assert.assertEquals("17.1", cgfCluster.getClusteredSpectra().get(0).getProperty(
                KnownProperties.ADDING_SCORE));
    }

    @Test
    public void testClusteringAppender() throws Exception {
        // get the spectra as clusters
        List<ICluster> clusters = ParserUtilities.readMGFClusters(testMgf);
        ICluster cluster = clusters.get(0);
        cluster.getClusteredSpectra().get(0).setProperty(KnownProperties.ADDING_SCORE, "17.1");

        StringBuilder stringBuilder = new StringBuilder();

        DotClusterClusterAppender.INSTANCE.appendCluster(stringBuilder, cluster);

        String clusteringResult = stringBuilder.toString();

        String[] lines = clusteringResult.split("\n");

        Assert.assertEquals(8, lines.length);

        String[] specFields = lines[7].split("\t");

        Assert.assertEquals(10, specFields.length);
        Assert.assertEquals("{\"RT\": \"2293.276122\", \"ADDING_SCORE\": \"17.1\"}", specFields[9]);
    }
}

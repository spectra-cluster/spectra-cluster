package uk.ac.ebi.pride.spectracluster.io;

import junit.framework.Assert;
import org.junit.Before;
import org.junit.Test;
import uk.ac.ebi.pride.spectracluster.cluster.ICluster;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.spectrum.KnownProperties;

import java.io.File;
import java.io.FileReader;
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
        StringBuilder stringBuilder = new StringBuilder();
        CGFClusterAppender.INSTANCE.appendCluster(stringBuilder, clusters.get(0));

        String testString = stringBuilder.toString();

        String[] lines = testString.split("\n");
        Assert.assertEquals(161, lines.length);
        Assert.assertEquals("RTINSECONDS=2293.276122", lines[8]);

        // read the cluster in again
        CGFSpectrumIterable clusterIterator = new CGFSpectrumIterable(new StringReader(testString));
        Assert.assertTrue(clusterIterator.iterator().hasNext());
        ICluster cgfCluster = clusterIterator.iterator().next();

        Assert.assertEquals(1, cgfCluster.getClusteredSpectraCount());
    }
}

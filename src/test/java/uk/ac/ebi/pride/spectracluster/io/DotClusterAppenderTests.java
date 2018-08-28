package uk.ac.ebi.pride.spectracluster.io;

import junit.framework.Assert;
import org.junit.Before;
import org.junit.Test;
import uk.ac.ebi.pride.spectracluster.cluster.GreedySpectralCluster;
import uk.ac.ebi.pride.spectracluster.cluster.ICluster;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.util.Defaults;

import java.io.File;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.List;

/**
 * Created by jg on 12.10.17.
 */
public class DotClusterAppenderTests {
    private List<ISpectrum> testSpectra;

    @Before
    public void setUp() throws Exception {
        Defaults.resetDefaults();
        File testFile = new File(DotClusterAppenderTests.class.getClassLoader().getResource("spectra_400.0_4.0.mgf").toURI());
        testSpectra = new ArrayList<ISpectrum>();
        ISpectrum[] readSpectra = ParserUtilities.readMGFScans(testFile);

        for (ISpectrum s : readSpectra)
            testSpectra.add(s);
    }

    @Test
    public void testCreateClusteringFile() throws Exception {
        StringWriter writer = new StringWriter(50 * 1024);
        DotClusterClusterAppender appender = new DotClusterClusterAppender(false);

        // create the cluster
        ICluster cluster = new GreedySpectralCluster("testId");

        for (int i = 0; i < 10; i++) {
            cluster.addSpectra(testSpectra.get(i));
        }

        // write teh cluster
        appender.appendStart(writer, "some version");
        appender.appendCluster(writer, cluster);
        appender.appendEnd(writer);

        // process the result as lines
        String[] resultLines = writer.toString().split("\n");

        Assert.assertEquals(25, resultLines.length);
        Assert.assertTrue(resultLines[12].startsWith("consensus_mz"));
        Assert.assertTrue(resultLines[13].startsWith("consensus_intens"));
        Assert.assertTrue(resultLines[14].startsWith("consensus_peak_counts"));

        // make sure all consensus lines have the same number of items
        int nMzValues = resultLines[12].split(",").length;
        int nIntensValues = resultLines[13].split(",").length;
        int nCountValues = resultLines[14].split(",").length;

        Assert.assertTrue(nMzValues == nIntensValues && nMzValues == nCountValues);

        // make sure the count values are not all 1
        String[] countValues = resultLines[14].split(",");
        boolean allOne = true;
        for (String stringValue : countValues) {
            if (stringValue != "1") {
                allOne = false;
                break;
            }
        }
        Assert.assertFalse(allOne);
    }
}
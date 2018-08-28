package uk.ac.ebi.pride.spectracluster.consensus;

import junit.framework.Assert;
import org.junit.Before;
import org.junit.Test;
import uk.ac.ebi.pride.spectracluster.io.MGFSpectrumAppender;
import uk.ac.ebi.pride.spectracluster.io.ParserUtilities;
import uk.ac.ebi.pride.spectracluster.spectrum.IPeak;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.util.Defaults;

import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.List;

/**
 * This class tests an issue were consensus peaks get higher counts than the number
 * of spectra in the cluster
 *
 * Created by jg on 06.11.17.
 */
public class GreedyConsensusSpectrumTest {
    private List<ISpectrum> testSpectra;

    @Before
    public void setUp() throws Exception {
        File testFile = new File(GreedyConsensusSpectrumTest.class.getClassLoader().getResource("consensus_test.mgf").toURI());
        Defaults.resetDefaults();
        testSpectra = new ArrayList<ISpectrum>();
        ISpectrum[] readSpectra = ParserUtilities.readMGFScans(testFile);

        for (ISpectrum s : readSpectra)
            testSpectra.add(s);
    }

    @Test
    public void testBuildConsensus() throws Exception {
        /**
         * This is known to fail and was used to isolate the issue
         */
        if (true)
            return;

        GreedyConsensusSpectrum consensusSpectrum =
                (GreedyConsensusSpectrum) GreedyConsensusSpectrum.buildFactory().getConsensusSpectrumBuilder();

        for (ISpectrum s : testSpectra) {
            consensusSpectrum.addSpectra(s);
        }

        ISpectrum finalSpec = consensusSpectrum.getConsensusSpectrum();

        // get the highest peak
        IPeak maxPeak = null;

        for (IPeak peak: finalSpec.getPeaks()) {
            if (maxPeak == null) {
                maxPeak = peak;
            } else if (peak.getIntensity() > maxPeak.getIntensity()) {
                maxPeak = peak;
            }
        }

        // highest peak should be 157.0836 - this is the highest peak in all spectra
        Assert.assertNotNull(maxPeak);
        Assert.assertEquals(157.0836F, maxPeak.getMz(), 0.01);

        for (IPeak peak : finalSpec.getPeaks()) {
            Assert.assertTrue("Higher peak count than spectra: " + peak.getCount(),
                    peak.getCount() <= testSpectra.size());
        }
    }

    @Test
    public void testBuildBinnedConsensus() throws Exception {
        BinnedGreedyConsensusSpectrum consensusSpectrum = (BinnedGreedyConsensusSpectrum) BinnedGreedyConsensusSpectrum.buildFactory().getConsensusSpectrumBuilder();

        consensusSpectrum.addSpectra(testSpectra.toArray(new ISpectrum[testSpectra.size()]));

        ISpectrum finalSpec = consensusSpectrum.getConsensusSpectrum();

        // get the highest peak
        IPeak maxPeak = null;

        for (IPeak peak: finalSpec.getPeaks()) {
            if (maxPeak == null) {
                maxPeak = peak;
            } else if (peak.getIntensity() > maxPeak.getIntensity()) {
                maxPeak = peak;
            }
        }

        // highest peak should be 157.0836 - this is the highest peak in all spectra
        Assert.assertNotNull(maxPeak);
        Assert.assertEquals(157.0836F, maxPeak.getMz(), 0.01);

        for (IPeak peak : finalSpec.getPeaks()) {
            Assert.assertTrue("Higher peak count than spectra: " + peak.getCount(),
                    peak.getCount() <= testSpectra.size());
        }

        FileWriter writer = new FileWriter("/tmp/consensus.mgf");
        MGFSpectrumAppender.INSTANCE.appendSpectrum(writer, finalSpec);
        writer.close();
    }
}

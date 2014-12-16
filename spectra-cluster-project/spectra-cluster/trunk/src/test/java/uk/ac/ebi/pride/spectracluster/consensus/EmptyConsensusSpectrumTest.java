package uk.ac.ebi.pride.spectracluster.consensus;

import junit.framework.Assert;
import org.junit.Test;
import uk.ac.ebi.pride.spectracluster.filter.IPeakFilter;
import uk.ac.ebi.pride.spectracluster.io.ParserUtilities;
import uk.ac.ebi.pride.spectracluster.spectrum.IPeak;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;

import java.io.File;

/**
 * Created by jg on 04.12.14.
 */
public class EmptyConsensusSpectrumTest {
    private final String testFilePath = "/home/jg/Projects/ebi-pride/pride-cluster-2/Testfiles/7d88f860-eec7-47b7-b574-0f9a6ef9ac09_mz.mgf";

    @Test
    public void testEmptySpectrum() {
        File inputFile = new File(testFilePath);

        // only run this test if the input file exists
        if (!inputFile.exists())
            return;

        ISpectrum[] mgfSpectra = ParserUtilities.readMGFScans(inputFile);

        IConsensusSpectrumBuilder consensusSpectrum = ConsensusSpectrum.buildFactory(IPeakFilter.NULL_FILTER).getConsensusSpectrumBuilder();

        consensusSpectrum.addSpectra(mgfSpectra);

        ISpectrum consensusSpec = consensusSpectrum.getConsensusSpectrum();

        for (IPeak peak : consensusSpec.getPeaks()) {
            Assert.assertTrue("Consensus spectrum contains peak with 0 m/z", peak.getMz() > 0);
            Assert.assertTrue("Consensus spectrum contains peak with 0 intensity", peak.getIntensity() > 0);
        }
    }
}

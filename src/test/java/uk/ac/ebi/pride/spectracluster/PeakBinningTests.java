package uk.ac.ebi.pride.spectracluster;

import org.junit.Test;
import uk.ac.ebi.pride.spectracluster.spectrum.IPeak;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.util.ClusteringTestUtilities;
import uk.ac.ebi.pride.spectracluster.util.ConsensusSpectraItems;
import uk.ac.ebi.pride.spectracluster.util.PeakUtilities;

import java.util.List;

/**
 * uk.ac.ebi.pride.spectracluster.PeakBinningTests
 *
 * @author Steve Lewis
 */
public class PeakBinningTests {


    @Test
    public void testPeakBinning() throws Exception {

        final List<ConsensusSpectraItems> consensusSpectraItems = ClusteringTestUtilities.readConsensusSpectraItemsFromResource();
        for (ConsensusSpectraItems csi : consensusSpectraItems) {
            for (ISpectrum sp : csi.getSpectra()) {
                testPeakBin(sp.getPeaks());
            }
        }

    }

    public void testPeakBin(List<IPeak> peaks) {
        int maxPerBin = 5;
        double minMZ = 0;
        double maxMZ = 5000;
        double binSize = 100;
        //noinspection UnusedDeclaration
        List<IPeak> binned = PeakUtilities.getHighestInBins(peaks, minMZ, maxMZ, binSize, maxPerBin);
//        for (IPeak iPeak : binned) {
//
//        }
    }
}

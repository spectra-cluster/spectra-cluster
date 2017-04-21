package uk.ac.ebi.pride.spectracluster.consensus;

import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;
import uk.ac.ebi.pride.spectracluster.spectrum.IPeak;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.util.*;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * @author Rui Wang
 * @version $Id$
 */
public class ConsensusSpectrumTests {

    private List<String> spectrumIds = new ArrayList<String>(Arrays.asList("83931", "1258781", "3722"));
    private List<ISpectrum> filteredOriginalSpectra = new ArrayList<ISpectrum>();


    @SuppressWarnings("UnusedDeclaration")
    public boolean peakListsEquivalent(List<IPeak> l1, List<IPeak> l2) {
        if (l1.size() != l2.size())
            return false;
        for (int i = 0; i < l1.size(); i++) {
            IPeak p1 = l1.get(i);
            IPeak p2 = l2.get(i);
            if (!p1.equivalent(p2))
                return false;
        }
        return true;
    }

    @Before
    public void setUp() throws Exception {
        List<ISpectrum> mgfSpectra = ClusteringTestUtilities.readISpectraFromResource();


        for (ISpectrum originalSpectrum : mgfSpectra) {
            if (spectrumIds.contains(originalSpectrum.getId())) {
                filteredOriginalSpectra.add(originalSpectrum);
            }
        }
    }


    @Test
    public void testConsensusSpectrum() throws Exception {

        IConsensusSpectrumBuilder consensusSpectrumBuilder = Defaults.getDefaultConsensusSpectrumBuilder();
        consensusSpectrumBuilder.addSpectra(filteredOriginalSpectra.toArray(new ISpectrum[filteredOriginalSpectra.size()]));

        ISpectrum consensusSpectrum1 = consensusSpectrumBuilder.getConsensusSpectrum();

        Collections.shuffle(filteredOriginalSpectra);
        consensusSpectrumBuilder =Defaults.getDefaultConsensusSpectrumBuilder();
        consensusSpectrumBuilder.addSpectra(filteredOriginalSpectra.toArray(new ISpectrum[filteredOriginalSpectra.size()]));
        ISpectrum consensusSpectrum2 = consensusSpectrumBuilder.getConsensusSpectrum();

        Collections.shuffle(filteredOriginalSpectra);
        consensusSpectrumBuilder = Defaults.getDefaultConsensusSpectrumBuilder();
        consensusSpectrumBuilder.addSpectra(filteredOriginalSpectra.toArray(new ISpectrum[filteredOriginalSpectra.size()]));
        ISpectrum consensusSpectrum3 = consensusSpectrumBuilder.getConsensusSpectrum();

        final boolean equivalent = consensusSpectrum1.equivalent(consensusSpectrum2);
        Assert.assertTrue(equivalent);
        final boolean equivalent1 = consensusSpectrum2.equivalent(consensusSpectrum3);
        Assert.assertTrue(equivalent1);
    }

//
//    @Test
//    public void testConsensusSpectrumStepByStep() throws Exception {
//        for (int i = 0; i < 10; i++) {
//            Collections.shuffle(filteredOriginalSpectra);
//            List<IPeak> allPeaks = consensusSpectrumBuilder.addAllPeaks(filteredOriginalSpectra);
//            List<IPeak> mergedIdenticalPeaks = consensusSpectrumBuilder.mergeIdenticalPeaks(allPeaks);
//
//
//            Collections.shuffle(filteredOriginalSpectra);
//            List<IPeak> allPeaks1 = consensusSpectrumBuilder.addAllPeaks(filteredOriginalSpectra);
//            Assert.assertTrue(peakListsEquivalent(allPeaks, allPeaks1));
//            List<IPeak> mergedIdenticalPeaks1 = consensusSpectrumBuilder.mergeIdenticalPeaks(allPeaks);
//            Assert.assertTrue(peakListsEquivalent(mergedIdenticalPeaks, mergedIdenticalPeaks1));
//        }
//    }
}

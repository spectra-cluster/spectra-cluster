package uk.ac.ebi.pride.spectracluster;

import org.junit.*;
import uk.ac.ebi.pride.spectracluster.cluster.*;
import uk.ac.ebi.pride.spectracluster.consensus.*;
import uk.ac.ebi.pride.spectracluster.filter.*;
import uk.ac.ebi.pride.spectracluster.io.*;
import uk.ac.ebi.pride.spectracluster.similarity.*;
import uk.ac.ebi.pride.spectracluster.spectrum.*;
import uk.ac.ebi.pride.spectracluster.util.*;

import java.io.*;
import java.util.*;

/**
 * uk.ac.ebi.pride.spectracluster.BinningPeakFilterTest
 *
 * @author Steve Lewis
 * @date 29/05/2014
 */
public class BinningPeakFilterTestMain {

    private static int numberOK = 0;
    private static int numberSlightlyOFF = 0;
    private static int numberOFF = 0;
    private static int numberWayOFF = 0;

    private static void processCGFFile(File file) {
        final ICluster[] clusters = ParserUtilities.readSpectralCluster(file);
        //noinspection ForLoopReplaceableByForEach
        for (int i = 0; i < clusters.length; i++) {
            ICluster cluster = clusters[i];
            testClusterConcensusSpectrum(cluster);
        }
    }

    private static void testClusterConcensusSpectrum(ICluster cluster) {
        final ISpectrum oldConcensusSpectrum = cluster.getConsensusSpectrum();
        final IConsensusSpectrumBuilder builder = ConsensusSpectrum.FACTORY.getConsensusSpectrumBuilder();
        for (ISpectrum spec : cluster.getClusteredSpectra()) {
            builder.addSpectra(spec);
        }
        ISpectrum newConcensusSpectrum = builder.getConsensusSpectrum();
        final ISimilarityChecker defaultSimilarityChecker = Defaults.getDefaultSimilarityChecker();
        double dotProduct = defaultSimilarityChecker.assessSimilarity(oldConcensusSpectrum, newConcensusSpectrum);
        if (dotProduct > 0.999)
            numberOK++;
        else {
            numberSlightlyOFF++;
            if (dotProduct < 0.9) {
                numberOFF++;
                validateCluster(cluster);
                if (dotProduct < 0.8) {
                    numberWayOFF++;
                    validateCluster(cluster);
                    System.out.println("Good " + numberOK + " off " + numberSlightlyOFF + " bad " + numberOFF + " very bad " + numberWayOFF);
                    //              Assert.assertTrue(numberWayOFF < 0.04 * numberOK);
                }

            }
        }


    }


    public static void validateCluster(ICluster cluster) {
        for (ISpectrum spectrum : cluster.getClusteredSpectra()) {
            final List<IPeak> oldPeaks = spectrum.getPeaks();
            final List<IPeak> newPeaks = BinnedHighestNPeakFilter.DEFAULT.filter(oldPeaks);
            testFilteredPeaks(oldPeaks, newPeaks);
        }
    }


    private static void testFilteredPeaks(List<IPeak> oldPeaks, List<IPeak> newPeaks) {
        for (int i = MZIntensityUtilities.LOWEST_USABLE_MZ; i < 5000; i += 50) {
            testPeaksInBin(i, oldPeaks, newPeaks);

        }

    }

    public static final int BIN_SIZE = BinnedHighestNPeakFilter.DEFAULT_BIN_SIZE;

    private static void testPeaksInBin(int bin, List<IPeak> oldPeaks, List<IPeak> newPeaks) {
        int oldCount = 0;
        int newCount = 0;
        double endBin = BIN_SIZE + bin;
        for (IPeak peak : oldPeaks) {
            double mz = peak.getMz();
            if (mz < bin)
                continue;
            if (mz > endBin)
                break;
            oldCount++;
        }
        for (IPeak peak : newPeaks) {
            double mz = peak.getMz();
            if (mz < bin)
                continue;
            if (mz > endBin)
                break;
            newCount++;
        }

        oldCount = Math.min(oldCount, BinnedHighestNPeakFilter.DEFAULT_MAX_PEAKS_PER_BIN);
        if (newCount < oldCount) {
            if (newCount < oldCount - 1) {
                if (newCount < oldCount - 2) {
                    if (newCount < oldCount - 3) {
                        Assert.assertTrue(newCount >= oldCount);
                    }
                }
            }
        }

    }


     public static void main(String[] args) {
        if (args.length == 0) {
            System.out.println("usage CGFFile CGFile ...");
        }
        //noinspection ForLoopReplaceableByForEach
        for (int i = 0; i < args.length; i++) {
            String arg = args[i];
            File udir = new File(arg);
            //noinspection ForLoopReplaceableByForEach
            for (File file : udir.listFiles()) {
                System.out.println(file);
                processCGFFile(file);
            }
        }

         MaximialPeakFilter.showStatistics(System.out);

         System.out.println("Good " + numberOK + " off " + numberSlightlyOFF + " bad " + numberOFF + " very bad " + numberWayOFF);
     }


}

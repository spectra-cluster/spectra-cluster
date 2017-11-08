package uk.ac.ebi.pride.spectracluster.consensus;

import uk.ac.ebi.pride.spectracluster.cluster.ICluster;
import uk.ac.ebi.pride.spectracluster.cluster.SpectralCluster;
import uk.ac.ebi.pride.spectracluster.io.ParserUtilities;
import uk.ac.ebi.pride.spectracluster.spectrum.IPeak;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.util.Defaults;

import java.io.File;

/**
 * uk.ac.ebi.pride.spectracluster.consensus.LargeConcensusSpectrumTest
 * User: Steve
 * Date: 9/11/2014
 */
public class LargeConcensusSpectrumTest {
    public static final LargeConcensusSpectrumTest[] EMPTY_ARRAY = {};

    /**
     * arg is a cgf with large spectra
     * @param args
     */
    public static void main(String[] args) {
        ICluster[] clusters = ParserUtilities.readSpectralCluster(new File(args[0]));

        SpectralCluster spectralCluster = new SpectralCluster("1234", Defaults.getDefaultConsensusSpectrumBuilder());
        /**
         * build ONE giant cluster
         */
        int totalSpectra = 0;
        long start = System.currentTimeMillis();
        for (int i = 0; i < clusters.length; i++) {
           long  addstart = System.currentTimeMillis();
            ICluster cluster = clusters[i];
            totalSpectra += cluster.getClusteredSpectraCount();
            for (ISpectrum sc : cluster.getClusteredSpectra()) {
                spectralCluster.addSpectra(sc);

            }
            long addEnd = System.currentTimeMillis();
            int elapsed = (int)((addEnd - addstart) / 1000);
            System.out.println("Total " + totalSpectra + " time " + elapsed);

        }
        ISpectrum consensusSpectrum = spectralCluster.getConsensusSpectrum();
        for (IPeak iPeak : consensusSpectrum.getPeaks()) {
            System.out.println(iPeak);
        }
        long totalEnd = System.currentTimeMillis();
         int elapsed = (int)((totalEnd - start) / 1000);
         System.out.println("Finished Total " + totalSpectra + " time " + elapsed);

    }
}

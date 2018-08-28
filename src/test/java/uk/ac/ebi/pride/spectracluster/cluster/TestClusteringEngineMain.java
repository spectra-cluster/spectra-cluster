package uk.ac.ebi.pride.spectracluster.cluster;


import uk.ac.ebi.pride.spectracluster.engine.IClusteringEngine;
import uk.ac.ebi.pride.spectracluster.io.CGFClusterAppender;
import uk.ac.ebi.pride.spectracluster.io.ParserUtilities;
import uk.ac.ebi.pride.spectracluster.util.Defaults;
import uk.ac.ebi.pride.tools.pride_spectra_clustering.impl.PrideClusteringEngine;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.List;

/**
 * uk.ac.ebi.pride.spectracluster.cluster.BinningClusteringMain
 * main to use BinningClusteringEngine - this is really a test main
 * User: Steve
 * Date: 7/5/13
 */
public class TestClusteringEngineMain {

    private long readTimeMillisec;

    public TestClusteringEngineMain() {
    }

    public long getReadTimeMillisec() {
        return readTimeMillisec;
    }

    public void addReadTimeMillisec(final long pReadTimeMillisec) {
        readTimeMillisec += pReadTimeMillisec;
    }

    /**
     * do the work of clustering one MGF file
     *
     * @param inputFile !null existing readable mgf file
     */
    protected void processFile(final File inputFile) {
        if (!inputFile.getName().toLowerCase().endsWith(".mgf"))
            return; // not an mgf

        long start = System.currentTimeMillis();
        List<ICluster> clusters = ParserUtilities.readMGFClusters(inputFile);

        /**
         * Add your favorite clustering engine here
         */
        //    PeakMatchClusteringEngine engine = new PeakMatchClusteringEngine();
        IClusteringEngine engine = new PrideClusteringEngine();
        long end = System.currentTimeMillis();
        final long readTIme = end - start;
        addReadTimeMillisec(readTIme);
        double seconds = (readTIme / 1000);
        //noinspection UnnecessaryLocalVariable,UnusedDeclaration,UnusedAssignment
        double min = seconds / 60;
        System.out.println("read " + inputFile + " with " + clusters.size() + " spectra in " + String.format("%10.3f", seconds).trim());

        if (clusters.size() == 0)     // nothing there
            return;
        if (clusters.size() == 1)    // no clustering to do
            return;
        start = System.currentTimeMillis();

        for (ICluster sc : clusters) {
            engine.addClusters(sc);
        }

        for (int i = 0; i < Defaults.getNumberReclusteringPasses(); i++) {
            if (!engine.processClusters()) {
                break;
            }
        }

        final List<ICluster> clusters1 = (List<ICluster>) engine.getClusters();


        saveClusters(clusters1, inputFile);


        end = System.currentTimeMillis();
        seconds = ((end - start) / 1000);
        min = seconds / 60;
        System.out.println("clustered " + inputFile + " with " + clusters.size() + " spectra in " + String.format("%10.3f sec", seconds).trim());
    }

    /**
     * write clusters to a file in the default directory with the extension .cgf
     *
     * @param pClusters1 !null list of clusters
     * @param pInputFile !null input file
     */
    protected void saveClusters(final List<ICluster> pClusters1, final File pInputFile) {

        if (pClusters1.size() == 0)
            return;
        String outName = pInputFile.getName().replace(".mgf", "") + ".cgf";
        PrintWriter out = null;
        try {
            System.out.println(new File(outName).getCanonicalPath());
            out = new PrintWriter(new FileWriter(outName));
            final CGFClusterAppender clusterAppender = CGFClusterAppender.INSTANCE;
            for (ICluster iPeptideSpectralCluster : pClusters1) {
                clusterAppender.appendCluster(out, iPeptideSpectralCluster);
            }
        } catch (IOException e) {
            throw new RuntimeException(e);

        } finally {
            if (out != null)
                out.close();
        }
    }

    /**
     * process every file in a directory containing mgf files
     *
     * @param pF !null existing directory
     */
    protected void processDirectory(final File pF) {
        long start = System.currentTimeMillis();
        long end = System.currentTimeMillis();
        double min;
        final File[] files = pF.listFiles();
        if (files != null) {
            for (File file : files) {
                if (file.isFile())
                    processFile(file); // todo better directory handling
                end = System.currentTimeMillis();
                int seconds = (int) ((end - start) / 1000);
                min = seconds / 60;
            }
        }
    }


    protected static void usage() {
        System.out.println("Usage <mgf file or directory> ...");
    }


    public static void main(String[] args) {
        if (args.length == 0) {
            usage();
            return;
        }


        TestClusteringEngineMain mainClusterer = new TestClusteringEngineMain();
        long start = System.currentTimeMillis();
        long end = System.currentTimeMillis();
        double min = 0;
        for (String arg : args) {
            File f = new File(arg);
            if (!f.exists())
                throw new IllegalArgumentException("File " + arg + " does not exist");
            if (f.isDirectory())
                mainClusterer.processDirectory(f);
            else
                mainClusterer.processFile(f);

            end = System.currentTimeMillis();
            int seconds = (int) ((end - start) / 1000);
            min = seconds / 60;
        }
        double readMin = mainClusterer.getReadTimeMillisec() / (60 * 1000);
        System.out.println("read in " + String.format("%10.2f", readMin) + " min");
        System.out.println("Processed in " + String.format("%10.2f", min - readMin) + " min");
        System.out.println("Total " + String.format("%10.2f", min) + " min");
    }


}

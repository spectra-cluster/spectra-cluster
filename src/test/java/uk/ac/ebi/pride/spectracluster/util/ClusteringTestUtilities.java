package uk.ac.ebi.pride.spectracluster.util;


import uk.ac.ebi.pride.spectracluster.cluster.FrankEtAClusterEngineTest;
import uk.ac.ebi.pride.spectracluster.cluster.ICluster;
import uk.ac.ebi.pride.spectracluster.cluster.SpectralCluster;
import uk.ac.ebi.pride.spectracluster.consensus.IConsensusSpectrumBuilder;
import uk.ac.ebi.pride.spectracluster.io.ParserUtilities;
import uk.ac.ebi.pride.spectracluster.similarity.ISimilarityChecker;
import uk.ac.ebi.pride.spectracluster.spectrum.IPeak;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.tools.pride_spectra_clustering.impl.Adapters;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.LineNumberReader;
import java.net.URISyntaxException;
import java.net.URL;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * uk.ac.ebi.pride.spectracluster.util.ClusteringTestUtilities
 * User: Steve
 * Date: 6/19/13
 */
public class ClusteringTestUtilities {

    public static final String SAMPLE_MGF_FILE = "spectra_400.0_4.0.mgf";
    public static final String SAMPLE_CLUSTERING_FILE = "Sample.clustering";

    public static final String SAMPLE_CGF_FILE = "uk/ac/ebi/pride/spectracluster/io/spectra_400.0_4.0.cgf";
    public static final String SAMPLE_SECOND_CGF_FILE = "uk/ac/ebi/pride/spectracluster/io/ClusterBin0400.cgf";

    /**
     * read a resource mgf as a list of spectra
     *
     * @return
     */
    public static List<ISpectrum> readConsensusSpectralItems() {
        return readConsensusSpectralItems(SAMPLE_MGF_FILE);
    }

    /**
     * read a resource mgf as a list of spectra
     *
     * @param resName
     * @return
     */
    public static List<ISpectrum> readConsensusSpectralItems(String resName) {
        File inputFile = getSpectrumFile(resName);


        ISpectrum[] mgfSpectra = ParserUtilities.readMGFScans(inputFile);
        return Arrays.asList(mgfSpectra);
    }

    public static File getSpectrumFile() {
        return getSpectrumFile(SAMPLE_MGF_FILE);
    }

    public static File getSpectrumFile(String resName) {
        // load a file contains a list of clusters
        URL url;
        url = ClusteringTestUtilities.class.getClassLoader().getResource(resName);
        if (url == null) {
            throw new IllegalStateException("no file for input found!");
        }
        File inputFile;
        try {
            inputFile = new File(url.toURI());
        } catch (URISyntaxException e) {
            throw new RuntimeException(e);

        }
        return inputFile;
    }

    //        // load a file contains a list of clusters
    //        URL url = FrankEtAClusterEngineTest.class.getClassLoader().getResource("uk/ac/ebi/pride/spectracluster/util/spectra_400.0_4.0.cgf");
    //        if (url == null) {
    //            throw new IllegalStateException("no file for input found!");
    //        }
    //        File inputFile = new File(url.toURI());
    //
    //        consensusSpectraItems = ParserUtilities.readClusters(inputFile);
    //
    //        url = FrankEtAClusterEngineTest.class.getClassLoader().getResource("uk/ac/ebi/pride/spectracluster/util/spectra_400.0_4.0.mgf");
    //        if (url == null) {
    //            throw new IllegalStateException("no file for input found!");
    //        }
    //        inputFile = new File(url.toURI());
    //
    //      //  consensusSpectraItems = ParserUtilities.guaranteeMGFParse(inputFile);

    /**
     * read a resource mgf as a list of spectra
     *
     * @return
     */
    public static List<ISpectrum> readISpectraFromResource() {
        return readISpectraFromResource(SAMPLE_MGF_FILE);
    }

    /**
     * read a resource mgf as a list of spectra
     *
     * @param resName
     * @return
     */
    public static List<ISpectrum> readISpectraFromResource(String resName) {
        // load a file contains a list of clusters
        File inputFile = getSpectrumFile(resName);

        ISpectrum[] mgfSpectra = ParserUtilities.readMGFScans(inputFile);
        return Arrays.asList(mgfSpectra);
    }


    public static LineNumberReader getResourceClusteringReader() {
        return getResourceClusteringReader(SAMPLE_CLUSTERING_FILE);
    }

    public static LineNumberReader getResourceClusteringReader(String resName) {
        URL testFile = ClusteringTestUtilities.class.getClassLoader().getResource(resName);

        assert testFile != null;
        try {
            File specFile = new File(testFile.toURI());
            return new LineNumberReader(new FileReader(specFile));
        } catch (URISyntaxException e) {
            throw new RuntimeException(e);
        } catch (FileNotFoundException e) {
            throw new RuntimeException(e);
        }

    }

    public static List<ICluster> readSpectraClustersFromResource() {
        return readSpectraClustersFromResource(SAMPLE_CGF_FILE);
    }


    public static List<ICluster> readSecondSpectraClustersFromResource() {
        return readSpectraClustersFromResource(SAMPLE_SECOND_CGF_FILE);
    }

    public static List<ICluster> readSpectraClustersFromResource(String resName) {
        List<ConsensusSpectraItems> items = readConsensusSpectraItemsFromResource(resName);
        int index = 1000;
        List<ICluster> holder = new ArrayList<ICluster>();
        for (ConsensusSpectraItems si : items) {
            ICluster cluster = new SpectralCluster(Integer.toString(index++), Defaults.getDefaultConsensusSpectrumBuilder());
            for (ISpectrum sr : si.getSpectra())
                cluster.addSpectra(sr);
            holder.add(cluster);
        }

        return holder;
    }

    /**
     * read a resource mgf as a list of spectra
     *
     * @return
     */
    public static List<ConsensusSpectraItems> readConsensusSpectraItemsFromResource() {
        return readConsensusSpectraItemsFromResource(SAMPLE_CGF_FILE);
    }

    /**
     * read a resource mgf as a list of spectra
     *
     * @param resName
     * @return
     */
    public static List<ConsensusSpectraItems> readConsensusSpectraItemsFromResource(String resName) {
        try {
            // load a file contains a list of clusters
            URL url = FrankEtAClusterEngineTest.class.getClassLoader().getResource(resName);
            if (url == null) {
                throw new IllegalStateException("no file for input found!");
            }
            File inputFile = new File(url.toURI());

            //noinspection UnusedDeclaration,UnnecessaryLocalVariable
            List<ConsensusSpectraItems> consensusSpectraItems = Arrays.asList(ParserUtilities.readClusters(inputFile));
            return consensusSpectraItems;

        } catch (URISyntaxException e) {
            throw new RuntimeException(e);

        }

    }

    @SuppressWarnings("UnusedDeclaration")
    public static final int MAX_PEAK_MISMATCHES = 3;
    public static final int NUMBER_HIGH_MATCHES = 5;
    public static final double FORGIVEN_MZ_DIFFERENCE = 0.2;
    public static final double FORGIVEN_INTENSITY_RATIO = 0.1;
    public static final double FORGIVEN_INTENSITY_DIFFERENCE = 200;

    /**
     * tolerate a few errors
     *
     * @param sp1
     * @param sp2
     * @return
     */
    public static boolean areSpectraVeryClose(ISpectrum sp1, ISpectrum sp2) {

        // easy case
        if (sp1.equivalent(sp2))
            return true;

        // Find mismatched peaks
        //noinspection MismatchedQueryAndUpdateOfCollection
        List<IPeak> mismatched = new ArrayList<IPeak>();
        double diff = Math.abs(sp1.getPrecursorMz() - sp2.getPrecursorMz());
        if (diff > MZIntensityUtilities.SMALL_MZ_DIFFERENCE)
            return false;

        if (sp1.getPeaksCount() != sp1.getPeaksCount())
            return false;

        int numberMisMatchedPeaks = 0;
        List<IPeak> l1 = sp1.getPeaks();
        List<IPeak> l2 = sp2.getPeaks();
        for (int i = 0; i < l1.size(); i++) {
            IPeak p1 = l1.get(i);
            IPeak p2 = l2.get(i);
            if (!p1.equivalent(p2)) {
                mismatched.add(p1);
                mismatched.add(p2);
                numberMisMatchedPeaks++;
            }
        }
        if (numberMisMatchedPeaks == 0)
            return true;  // probably will not happen

        final ISimilarityChecker checker = Defaults.getDefaultSimilarityChecker();


        // well we better agree on the highest peaks
        final ISpectrum hp1 = sp1.getHighestNPeaks(NUMBER_HIGH_MATCHES);
        final ISpectrum hp2 = sp2.getHighestNPeaks(NUMBER_HIGH_MATCHES);
        l1 = hp1.getPeaks();
        l2 = hp2.getPeaks();
        for (int i = 0; i < l1.size(); i++) {
            IPeak p1 = l1.get(i);
            IPeak p2 = l2.get(i);
            double diffx = Math.abs(p1.getMz() - p2.getMz());
            // allow larger MZ diffrences
            if (diffx > FORGIVEN_MZ_DIFFERENCE)
                return false;

            double avgIntensity1 = p1.getIntensity() / p1.getCount();
            double avgIntensity2 = p2.getIntensity() / p2.getCount();
            double totalIntensity = avgIntensity1 + avgIntensity2;
            double diffI = Math.abs((avgIntensity1 / totalIntensity) - (avgIntensity2 / totalIntensity));
            // allow larger MZ diffrences
            if (diffI > FORGIVEN_INTENSITY_RATIO) {
                if (Math.abs(avgIntensity1 - avgIntensity2) > FORGIVEN_INTENSITY_DIFFERENCE)
                    return false;
            }

        }
        // ignore count and intensity diferences
        // dot product
        double similarity = checker.assessSimilarity(sp1, sp2);
        // certainly in the same cluster
        double defaultThreshold = Defaults.getSimilarityThreshold();  // 0.6 - 0.7
        //noinspection UnusedDeclaration
        double highThreshold = 1.0 - ((1.0 - defaultThreshold) / 2); // 0.8 - 0.85
        //noinspection RedundantIfStatement
        if (similarity < defaultThreshold)     // similarity better be 0.95
            return false;  // dot producta not close enough

        return true;

    }


    @SuppressWarnings("UnusedDeclaration")
    public static boolean peakListsEquivalent(List<IPeak> l1, List<IPeak> l2) {
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

    @SuppressWarnings("UnusedDeclaration")
    public static boolean mzEquivalent(List<IPeak> l1, List<uk.ac.ebi.pride.tools.pride_spectra_clustering.util.Peak> l2) {
        if (l1.size() != l2.size())
            return false;

        List<Float> mzValues1 = new ArrayList<Float>(l1.size());
        for (IPeak p : l1)
            mzValues1.add((float) MZIntensityUtilities.round(p.getMz(), 10)); // don't bother about rounding differences

        List<Float> mzValues2 = new ArrayList<Float>(l2.size());
        for (uk.ac.ebi.pride.tools.pride_spectra_clustering.util.Peak p : l2)
            mzValues2.add((float) MZIntensityUtilities.round(p.getMz(), 10));

        Collections.sort(mzValues1);
        Collections.sort(mzValues2);

        // print out all m/z values
//        for (int i = 0; i < mzValues1.size(); i++) {
//            System.out.println("Old " + i + "\t" + mzValues2.get(i) + "\tNem " + i + "\t" + mzValues1.get(i));
//        }

        for (int i = 0; i < mzValues1.size(); i++) {
            float mz1 = mzValues1.get(i);
            float mz2 = mzValues2.get(i);

            if (mz1 != mz2)
                return false;
        }

        return true;
    }

    public static boolean arePeakListsEquivalent(List<IPeak> l1x, List<uk.ac.ebi.pride.tools.pride_spectra_clustering.util.Peak> l2x) {

        int l1Count = PeakUtilities.getTotalCount(l1x);   // for debugging
        int l2Count = Adapters.getTotalCount(l2x);   // for debugging


        boolean isEqual = true;

        if (l1x.size() != l2x.size())
            return false;

        // copy and sort by
        List<IPeak> l1 = new ArrayList<IPeak>(l1x);
        Collections.sort(l1);
        List<IPeak> l2 = Adapters.fromPeaks(l2x);
        Collections.sort(l2);

        for (int i = 0; i < l1.size(); i++) {
            IPeak p1 = l1.get(i);
            IPeak p2 = l2.get(i);

            final int count1 = p1.getCount();
            final int count2 = p2.getCount();
            System.out.format(i + ": new = old\tm/z: %f = %f\t\tintens: %f = %f\tcount: %d = %d", p1.getMz(), p2.getMz(), p1.getIntensity(), p2.getIntensity(), count1, count2);

            if (count1 != count2) {
                System.out.println(" <-- count differs!");
                isEqual = false;
            } else if (MZIntensityUtilities.round(p1.getMz()) != MZIntensityUtilities.round(p2.getMz())) {
                System.out.println(" <-- m/z differ!");
                isEqual = false;
            } else if (MZIntensityUtilities.round(p1.getIntensity(), 100) != MZIntensityUtilities.round(p2.getIntensity(), 100)) {
                if (Math.abs(p1.getIntensity() - p2.getIntensity()) >= 2 || p1.getIntensity() < 100000) {
                    System.out.println(" <-- intensity differ!");
                    isEqual = false;
                } else {
                    System.out.println("");
                }
            } else {
                System.out.println("");
            }
        }
        // not same total
        //noinspection SimplifiableIfStatement
        if (l1Count != l2Count)
            return false;
        return isEqual;
    }

    public static boolean areNewPeakListsEquivalent(List<IPeak> l1, List<IPeak> l2, boolean print) {
        boolean isEqual = true;

        for (int i = 0; i < l1.size(); i++) {
            IPeak p1 = l1.get(i);
            IPeak p2 = l2.get(i);

            if (print)
                System.out.format(i + ": new = old\tm/z: %f = %f\t\tintens: %f = %f\tcount: %d = %d", p1.getMz(), p2.getMz(), p1.getIntensity(), p2.getIntensity(), p1.getCount(), p2.getCount());

            if (p1.getCount() != p2.getCount()) {
                if (print) System.out.println(" <-- count differs!");
                isEqual = false;
            } else if (p1.getMz() != p2.getMz()) {
                if (print) System.out.println(" <-- m/z differ!");
                isEqual = false;
            } else if (p1.getIntensity() != p2.getIntensity()) {
                if (print) System.out.println(" <-- intensity differ!");
                isEqual = false;
            } else {
                if (print) System.out.println("");
            }
        }
        return isEqual;
    }

    /**
     * create a list of consensusSpectra from a list of clusters
     *
     * @param pClusters !null cluster list
     * @return !null list of  consensusSpectra
     */
    public static List<ISpectrum> buildConsessusSpectra(final List<ICluster> pClusters, final IConsensusSpectrumBuilder consensusSpectrumBuilder) {
        List<ISpectrum> holder = new ArrayList<ISpectrum>();
        for (ICluster cluster : pClusters) {

            final List<ISpectrum> css = cluster.getClusteredSpectra();
            for (ISpectrum cs : css) {
                consensusSpectrumBuilder.addSpectra(cs);
            }

            final ISpectrum oldSpec = consensusSpectrumBuilder.getConsensusSpectrum();
            holder.add(oldSpec);
        }


        return holder;
    }


}



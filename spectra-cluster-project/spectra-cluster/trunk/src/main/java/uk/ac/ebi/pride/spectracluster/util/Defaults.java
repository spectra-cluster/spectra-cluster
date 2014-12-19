package uk.ac.ebi.pride.spectracluster.util;

import uk.ac.ebi.pride.spectracluster.consensus.*;
import uk.ac.ebi.pride.spectracluster.engine.EngineFactories;
import uk.ac.ebi.pride.spectracluster.engine.IClusteringEngine;
import uk.ac.ebi.pride.spectracluster.filter.BinnedHighestNPeakFilter;
import uk.ac.ebi.pride.spectracluster.filter.IPeakFilter;
import uk.ac.ebi.pride.spectracluster.filter.MaximialPeakFilter;
import uk.ac.ebi.pride.spectracluster.normalizer.IIntensityNormalizer;
import uk.ac.ebi.pride.spectracluster.normalizer.TotalIntensityNormalizer;
import uk.ac.ebi.pride.spectracluster.quality.IQualityScorer;
import uk.ac.ebi.pride.spectracluster.quality.SignalToNoiseChecker;
import uk.ac.ebi.pride.spectracluster.similarity.FrankEtAlDotProduct;
import uk.ac.ebi.pride.spectracluster.similarity.ISimilarityChecker;
import uk.ac.ebi.pride.spectracluster.util.comparator.ClusterComparator;

/**
 * uk.ac.ebi.pride.spectracluster.util.Defaults
 *
 * @author Steve Lewis
 * @date 20/05/13
 */
public class Defaults {

    public static final int DEFAULT_NUMBER_RECLUSTERING_PASSES = 2;

    public static final int DEFAULT_NUMBER_COMPARED_PEAKS = 15;

    public static final double DEFAULT_MZ_RANGE = 0.5;

    /**
     * Defines the similarity threshold above which spectra are added
     * to a cluster.
     */
    public static final double DEFAULT_SIMILARITY_THRESHOLD = 0.7;

    /**
     * Defines the similarity threshold below which spectra are removed
     * from a cluster.
     */
    public static final double DEFAULT_RETAIN_THRESHOLD = 0.6;

    public static final int DEFAULT_LARGE_BINNING_REGION = 1000;

    private static double similarityThreshold = DEFAULT_SIMILARITY_THRESHOLD;

    private static int largeBinningRegion = DEFAULT_LARGE_BINNING_REGION;

    private static int numberComparedPeaks = DEFAULT_NUMBER_COMPARED_PEAKS;

    private static double similarityMZRange = DEFAULT_MZ_RANGE;

    private static double retainThreshold = DEFAULT_RETAIN_THRESHOLD;

    private static int numberReclusteringPasses = DEFAULT_NUMBER_RECLUSTERING_PASSES;

    private static final int DEFAULT_MAJOR_PEAKS = 6;

    private static int majorPeakCount = DEFAULT_MAJOR_PEAKS;

    public static int getMajorPeakCount() {
        return majorPeakCount;
    }

    public static void setMajorPeakCount(int majorPeakCount) {
        Defaults.majorPeakCount = majorPeakCount;
    }

    public static double getSimilarityThreshold() {
        return similarityThreshold;
    }

    public static int getLargeBinningRegion() {
        return largeBinningRegion;
    }

    public static int getNumberComparedPeaks() {
        return numberComparedPeaks;
    }

    public static double getSimilarityMZRange() {
        return similarityMZRange;
    }

    public static double getRetainThreshold() {
        return retainThreshold;
    }

    public static void setSimilarityThreshold(double similarityThreshold) {
        Defaults.similarityThreshold = similarityThreshold;
    }

    public static void setLargeBinningRegion(int largeBinningRegion) {
        Defaults.largeBinningRegion = largeBinningRegion;
    }

    public static void setNumberComparedPeaks(int numberComparedPeaks) {
        Defaults.numberComparedPeaks = numberComparedPeaks;
    }

    public static void setSimilarityMZRange(double similarityMZRange) {
        Defaults.similarityMZRange = similarityMZRange;
    }

    /**
     * The retain threshold defines the similarity threshold below which
     * spectra are removed from a cluster.
     * @param retainThreshold
     */
    public static void setRetainThreshold(double retainThreshold) {
        Defaults.retainThreshold = retainThreshold;
    }

    public static int getNumberReclusteringPasses() {
        return numberReclusteringPasses;
    }

    public static void setNumberReclusteringPasses(final int pNumberReclusteringPasses) {
        numberReclusteringPasses = pNumberReclusteringPasses;
    }

    private static IConsensusSpectrumBuilder defaultConsensusSpectrumBuilder = null;

    /**
      * filter sees that we dont pass more then MaximialPeakFilter.DEFAULT_MAX_PEAKS peaks (100)
      */
   //  public static IPeakFilter defaultPeakFilter = IPeakFilter.NULL_FILTER; //   Take out for testing onlyu
     // public static IPeakFilter defaultPeakFilter = new MaximialPeakFilter(MaximialPeakFilter.DEFAULT_MAX_PEAKS); jg - this setting was active until 16-Dec-2014
    public static IPeakFilter defaultPeakFilter = new BinnedHighestNPeakFilter(20, 100, 50); // keep 20 peaks per 100 m/z with a 50 m/z overlap

     public static IPeakFilter getDefaultPeakFilter() {
         return defaultPeakFilter;
     }

     public static void setDefaultPeakFilter(IPeakFilter defaultPeakFilter) {
         Defaults.defaultPeakFilter = defaultPeakFilter;
     }


    /**
     * filter to use a consensus spectrum
     */
    private static ConcensusSpectrumBuilderFactory consensusFactory =  ConsensusSpectrum.FACTORY;

    public static ConcensusSpectrumBuilderFactory getConsensusFactory() {
        return consensusFactory;
    }

    public static void setConsensusFactory(ConcensusSpectrumBuilderFactory consensusFactory) {
        Defaults.consensusFactory = consensusFactory;
    }

    /**
     * this is the way to get a ConsensusSpectrumBuilder
     *
     * @return
     */
    public static IConsensusSpectrumBuilder getDefaultConsensusSpectrumBuilder() {
        return consensusFactory.getConsensusSpectrumBuilder();
    }


    private static ISimilarityChecker defaultSimilarityChecker = new FrankEtAlDotProduct(getSimilarityMZRange(), getNumberComparedPeaks());

    public static ISimilarityChecker getDefaultSimilarityChecker() {
        return defaultSimilarityChecker;
    }

    public static void setDefaultSimilarityChecker(ISimilarityChecker defaultSimilarityChecker) {
        Defaults.defaultSimilarityChecker = defaultSimilarityChecker;
    }

    private static IQualityScorer defaultQualityScorer = new SignalToNoiseChecker();

    public static IQualityScorer getDefaultQualityScorer() {
        return defaultQualityScorer;
    }

    public static void setDefaultQualityScorer(IQualityScorer defaultQualityScorer) {
        Defaults.defaultQualityScorer = defaultQualityScorer;
    }

    private static ClusterComparator defaultSpectrumComparator = ClusterComparator.INSTANCE;

    public static ClusterComparator getDefaultSpectrumComparator() {
        return defaultSpectrumComparator;
    }

    public static void setDefaultSpectrumComparator(ClusterComparator dc) {
        defaultSpectrumComparator = dc;
    }

     public static IClusteringEngine getDefaultClusteringEngine() {
        ISimilarityChecker similarityChecker = getDefaultSimilarityChecker();
        final ClusterComparator spectrumComparator = getDefaultSpectrumComparator();
        final double st = getSimilarityThreshold();
        final double rt = getRetainThreshold();
        return EngineFactories.buildClusteringEngineFactory(similarityChecker, spectrumComparator, st, rt).buildInstance();

    }

    /**
     * Default intensity normalizer
     */
    private static IIntensityNormalizer defaultIntensityNormalizer = TotalIntensityNormalizer.DEFAULT;

    public static IIntensityNormalizer getDefaultIntensityNormalizer() {
        return defaultIntensityNormalizer;
    }

    public static void setDefaultIntensityNormalizer(IIntensityNormalizer defaultIntensityNormalizer) {
        Defaults.defaultIntensityNormalizer = defaultIntensityNormalizer;
    }
}

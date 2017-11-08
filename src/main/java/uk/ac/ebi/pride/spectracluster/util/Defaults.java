package uk.ac.ebi.pride.spectracluster.util;

import uk.ac.ebi.pride.spectracluster.cdf.CumulativeDistributionFunction;
import uk.ac.ebi.pride.spectracluster.cdf.INumberOfComparisonAssessor;
import uk.ac.ebi.pride.spectracluster.cdf.MinNumberComparisonsAssessor;
import uk.ac.ebi.pride.spectracluster.consensus.ConcensusSpectrumBuilderFactory;
import uk.ac.ebi.pride.spectracluster.consensus.ConsensusSpectrum;
import uk.ac.ebi.pride.spectracluster.consensus.IConsensusSpectrumBuilder;
import uk.ac.ebi.pride.spectracluster.engine.EngineFactories;
import uk.ac.ebi.pride.spectracluster.engine.IClusteringEngine;
import uk.ac.ebi.pride.spectracluster.normalizer.IIntensityNormalizer;
import uk.ac.ebi.pride.spectracluster.normalizer.TotalIntensityNormalizer;
import uk.ac.ebi.pride.spectracluster.quality.IQualityScorer;
import uk.ac.ebi.pride.spectracluster.quality.SignalToNoiseChecker;
import uk.ac.ebi.pride.spectracluster.similarity.CombinedFisherIntensityTest;
import uk.ac.ebi.pride.spectracluster.similarity.ISimilarityChecker;
import uk.ac.ebi.pride.spectracluster.spectrum.IPeak;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.util.comparator.ClusterComparator;
import uk.ac.ebi.pride.spectracluster.util.function.Functions;
import uk.ac.ebi.pride.spectracluster.util.function.IFunction;
import uk.ac.ebi.pride.spectracluster.util.function.peak.FractionTICPeakFunction;
import uk.ac.ebi.pride.spectracluster.util.function.spectrum.HighestNSpectrumPeaksFunction;
import uk.ac.ebi.pride.spectracluster.util.function.spectrum.RemoveImpossiblyHighPeaksFunction;
import uk.ac.ebi.pride.spectracluster.util.function.spectrum.RemovePrecursorPeaksFunction;

import java.util.List;

/**
 * uk.ac.ebi.pride.spectracluster.util.Defaults
 *
 * @author Steve Lewis
 */
public class Defaults {

    public static final int DEFAULT_NUMBER_RECLUSTERING_PASSES = 2;

    public static final int DEFAULT_NUMBER_COMPARED_PEAKS = 15;

    public static final float DEFAULT_FRAGMENT_ION_TOLERANCE = 0.5F;

    /**
     * This default precursor tolerance is used by the incremental
     * clustering engines as window size
     */
    public static final float DEFAULT_PRECURSOR_ION_TOLERANCE = 2.0F;

    /**
     * Defines the similarity threshold above which spectra are added
     * to a cluster.
     */
    public static final double DEFAULT_SIMILARITY_THRESHOLD = 0.99;
    //public static final double DEFAULT_SIMILARITY_THRESHOLD = 48;

    /**
     * Defines the similarity threshold below which spectra are removed
     * from a cluster.
     * This is not used in greedy clustering
     */
    public static final double DEFAULT_RETAIN_THRESHOLD = 0.6; // not used in greedy clustering
    //public static final double DEFAULT_RETAIN_THRESHOLD = 45;

    public static final int DEFAULT_LARGE_BINNING_REGION = 1000;

    /**
     * If this is set, this minimum number of comparisons to calculate
     * the probabilities of a match based on the corresponding CDF
     * <p>
     * This is currently only implemented by the GreedyIncrementalClusteringEngine.
     */
    public static final INumberOfComparisonAssessor DEFAULT_NUMBER_COMPARISON_ASSESSOR = new MinNumberComparisonsAssessor(0);
    /**
     * By default no default cumulative distribution function is set. In this
     * case the CDF is loaded from the matching resources (see
     * CumulativeDistributionFunctionFactory)
     */
    public static final CumulativeDistributionFunction DEFAULT_CUMULATIVE_DISTRIBUTION_FUNCTION = null;

    private static double similarityThreshold = DEFAULT_SIMILARITY_THRESHOLD;

    private static int largeBinningRegion = DEFAULT_LARGE_BINNING_REGION;

    private static int numberComparedPeaks = DEFAULT_NUMBER_COMPARED_PEAKS;

    private static float fragmentIonTolerance = DEFAULT_FRAGMENT_ION_TOLERANCE;

    private static double retainThreshold = DEFAULT_RETAIN_THRESHOLD;

    private static int numberReclusteringPasses = DEFAULT_NUMBER_RECLUSTERING_PASSES;

    private static final int DEFAULT_MAJOR_PEAKS = 5;

    private static int majorPeakCount = DEFAULT_MAJOR_PEAKS;

    private static float defaultPrecursorIonTolerance = DEFAULT_PRECURSOR_ION_TOLERANCE;

    private static INumberOfComparisonAssessor numberOfComparisonAssessor = DEFAULT_NUMBER_COMPARISON_ASSESSOR;

    private static CumulativeDistributionFunction cumulativeDistributionFunction = DEFAULT_CUMULATIVE_DISTRIBUTION_FUNCTION;

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

    public static float getFragmentIonTolerance() {
        return fragmentIonTolerance;
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

    public static void setFragmentIonTolerance(float fragmentIonTolerance) {
        Defaults.fragmentIonTolerance = fragmentIonTolerance;

        Defaults.updateFragmentToleranceDependentDefaults();
    }

    /**
     * Update all default parameters that depend on the
     * fragment ion tolerance.
     *
     * Parameters that were set by the user are not updated
     *
     * In case of the DefaultSimilarityChecker, the
     * setFragmentIonTolerance is called. Therefore, this threshold
     * is automatically adapted.
     */
    private static void updateFragmentToleranceDependentDefaults() {
        if (!Defaults.defaultPeakFilterIsOverwritten) {
            Defaults.defaultPeakFilter = Defaults.generateDefaultPeakFilter();
        }

        Defaults.defaultSimilarityChecker.setFragmentIonTolerance(Defaults.fragmentIonTolerance);
    }

    public static float getDefaultPrecursorIonTolerance() {
        return defaultPrecursorIonTolerance;
    }

    public static void setDefaultPrecursorIonTolerance(float defaultPrecursorIonTolerance) {
        Defaults.defaultPrecursorIonTolerance = defaultPrecursorIonTolerance;
    }

    public static INumberOfComparisonAssessor getNumberOfComparisonAssessor() {
        return numberOfComparisonAssessor;
    }

    public static void setNumberOfComparisonAssessor(INumberOfComparisonAssessor numberOfComparisonAssessor) {
        Defaults.numberOfComparisonAssessor = numberOfComparisonAssessor;
    }

    /**
     * The retain threshold defines the similarity threshold below which
     * spectra are removed from a cluster.
     *
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

    /**
     * Generates the defaultPeakFiltering filter which is applied to spectra
     * right after they are loaded from the peak list file.
     *
     * Some of these filters depend on the set fragmentIonTolerance. Therefore,
     * they have to be regenerated when the tolerance is changed.
     * @return The default peak filtering function
     */
    private static IFunction<ISpectrum, ISpectrum> generateDefaultPeakFilter() {
        return Functions.join(new RemoveImpossiblyHighPeaksFunction(),
                Functions.join(
                        new RemovePrecursorPeaksFunction(fragmentIonTolerance),
                        new HighestNSpectrumPeaksFunction(150)));
    }

    /**
     * Indicates whether the user overwrote the default peak filtering
     * function. In this case, it should not be regenerated if the
     * fragment tolerance is changed.
     */
    private static boolean defaultPeakFilterIsOverwritten = false;

    private static IFunction<ISpectrum, ISpectrum> defaultPeakFilter = generateDefaultPeakFilter();

    public static IFunction<ISpectrum, ISpectrum> getDefaultPeakFilter() {
        return defaultPeakFilter;
    }

    public static void setDefaultPeakFilter(IFunction<ISpectrum, ISpectrum> defaultPeakFilter) {
        Defaults.defaultPeakFilter = defaultPeakFilter;
        Defaults.defaultPeakFilterIsOverwritten = true;
    }

    /**
     * default filter to use before comparing two spectra
     */
    private static IFunction<List<IPeak>, List<IPeak>> defaultComparisonPeakFilter =
            new FractionTICPeakFunction(0.5F, 20);

    public static IFunction<List<IPeak>, List<IPeak>> getDefaultComparisonPeakFilter() {
        return defaultComparisonPeakFilter;
    }

    public static void setDefaultComparisonPeakFilter(IFunction<List<IPeak>, List<IPeak>> defaultComparisonPeakFilter) {
        Defaults.defaultComparisonPeakFilter = defaultComparisonPeakFilter;
    }

    /**
     * filter to use a consensus spectrum
     */
    private static ConcensusSpectrumBuilderFactory consensusFactory = ConsensusSpectrum.FACTORY;

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

    //private static ISimilarityChecker defaultSimilarityChecker = new FrankEtAlDotProduct(getFragmentIonTolerance(), getNumberComparedPeaks());
    //private static ISimilarityChecker defaultSimilarityChecker = new FisherExactTest((float) getFragmentIonTolerance());
    private static ISimilarityChecker defaultSimilarityChecker = new CombinedFisherIntensityTest((float) getFragmentIonTolerance());

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

    /**
     * The default cumulative distribution function. If this is not
     * set, NULL is returned. In this case the matching CDF should be
     * fetched using the CumulativeDistributionFunctionFactory function
     * to load the CDF from the matching resource file.
     *
     * @return
     */
    public static CumulativeDistributionFunction getCumulativeDistributionFunction() {
        return cumulativeDistributionFunction;
    }

    /**
     * Sets the cumulative distribution function to be used. This overwrites
     * the cumulative distirbution function which is otherwise loaded from the
     * matching resources.
     *
     * @param cumulativeDistributionFunction
     */
    public static void setCumulativeDistributionFunction(CumulativeDistributionFunction cumulativeDistributionFunction) {
        Defaults.cumulativeDistributionFunction = cumulativeDistributionFunction;
    }

    /**
     * The minimum number of peaks in a consensus spectrum before the
     * peak filtering is used (retaining N peaks per M m/z).
     */
    public static final int DEFAULT_CONSENSUS_MIN_PEAKS = 50;

    private static int defaultConsensusMinPeaks = DEFAULT_CONSENSUS_MIN_PEAKS;

    /**
     * The minimum number of peaks in a consensus spectrum before the
     * peak filtering is used (retaining N peaks per M m/z).
     *
     * @return
     */
    public static int getDefaultConsensusMinPeaks() {
        return defaultConsensusMinPeaks;
    }

    /**
     * The minimum number of peaks in a consensus spectrum before the
     * peak filtering is used (retaining N peaks per M m/z).
     * @param defaultConsensusMinPeaks
     */
    public static void setDefaultConsensusMinPeaks(int defaultConsensusMinPeaks) {
        Defaults.defaultConsensusMinPeaks = defaultConsensusMinPeaks;
    }

    /**
     * Indicates whether additional debug information is stored in the spectra
     * objects during clustering (f.e. the number of comparisosn when adding a
     * spectrum).
     */
    private static boolean saveDebugInformation = false;

    /**
     * Indicates whether additional debug information should be saved
     * @return
     */
    public static boolean isSaveDebugInformation() {
        return saveDebugInformation;
    }

    /**
     * Sets whether additional debug information should be saved during clustering
     * @param saveDebugInformation
     */
    public static void setSaveDebugInformation(boolean saveDebugInformation) {
        Defaults.saveDebugInformation = saveDebugInformation;
    }

    /**
     * If set, the similarity score that led to the spectrum being added to
     * the cluster is saved as a property of each spectrum
     */
    private static boolean saveAddingScore = false;

    /**
     * Indicates whether the similarity score a spectrum has when it
     * was added to a cluster is saved.
     * @return
     */
    public static boolean isSaveAddingScore() {
        return saveAddingScore;
    }

    /**
     * Set whether the similarity score a spectrum has when being
     * added to a cluster should be saved.
     * @param saveAddingScore
     */
    public static void setSaveAddingScore(boolean saveAddingScore) {
        Defaults.saveAddingScore = saveAddingScore;
    }

    /**
     * Reset all values to their defaults
     */
    public static void resetDefaults() {
        defaultComparisonPeakFilter = new FractionTICPeakFunction(0.5F, 20);
        defaultConsensusMinPeaks = DEFAULT_CONSENSUS_MIN_PEAKS;
        defaultIntensityNormalizer = TotalIntensityNormalizer.DEFAULT;
        defaultPeakFilter = generateDefaultPeakFilter();
        defaultPeakFilterIsOverwritten = false;
        defaultPrecursorIonTolerance = DEFAULT_PRECURSOR_ION_TOLERANCE;
        defaultQualityScorer = new SignalToNoiseChecker();
        defaultSimilarityChecker = new CombinedFisherIntensityTest((float) getFragmentIonTolerance());
        defaultSpectrumComparator = ClusterComparator.INSTANCE;
        fragmentIonTolerance = DEFAULT_FRAGMENT_ION_TOLERANCE;
        largeBinningRegion = DEFAULT_LARGE_BINNING_REGION;
        majorPeakCount = DEFAULT_MAJOR_PEAKS;
        numberOfComparisonAssessor = DEFAULT_NUMBER_COMPARISON_ASSESSOR;
        numberComparedPeaks = DEFAULT_NUMBER_COMPARED_PEAKS;
        numberReclusteringPasses = DEFAULT_NUMBER_RECLUSTERING_PASSES;
        retainThreshold = DEFAULT_RETAIN_THRESHOLD;
        saveAddingScore = false;
        saveDebugInformation = false;
        similarityThreshold = DEFAULT_SIMILARITY_THRESHOLD;
    }
}

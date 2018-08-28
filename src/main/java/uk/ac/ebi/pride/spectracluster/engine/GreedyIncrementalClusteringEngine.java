package uk.ac.ebi.pride.spectracluster.engine;

import uk.ac.ebi.pride.spectracluster.cdf.CumulativeDistributionFunction;
import uk.ac.ebi.pride.spectracluster.cdf.CumulativeDistributionFunctionFactory;
import uk.ac.ebi.pride.spectracluster.cdf.INumberOfComparisonAssessor;
import uk.ac.ebi.pride.spectracluster.cluster.GreedySpectralCluster;
import uk.ac.ebi.pride.spectracluster.cluster.ICluster;
import uk.ac.ebi.pride.spectracluster.similarity.ISimilarityChecker;
import uk.ac.ebi.pride.spectracluster.spectrum.IPeak;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.spectrum.KnownProperties;
import uk.ac.ebi.pride.spectracluster.spectrum.Spectrum;
import uk.ac.ebi.pride.spectracluster.util.Defaults;
import uk.ac.ebi.pride.spectracluster.util.MZIntensityUtilities;
import uk.ac.ebi.pride.spectracluster.util.NumberUtilities;
import uk.ac.ebi.pride.spectracluster.util.function.IFunction;
import uk.ac.ebi.pride.spectracluster.util.predicate.IComparisonPredicate;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

/**
 * uk.ac.ebi.pride.spectracluster.engine.IncrementalClusteringEngine
 * a version of a clustering engine in which spectra are added incrementally and
 * clusters are shed when they are too far to use
 * User: Steve
 * Date: 7/5/13
 */
public class GreedyIncrementalClusteringEngine implements IIncrementalClusteringEngine {
    private final List<GreedySpectralCluster> clusters = new ArrayList<>();
    private final List<ISpectrum> filteredConsensusSpectra = new ArrayList<>();

    private final ISimilarityChecker similarityChecker;
    private final Comparator<ICluster> spectrumComparator;
    private final double windowSize;
    private final double mixtureProbability;
    private final CumulativeDistributionFunction cumulativeDistributionFunction;
    private final IFunction<List<IPeak>, List<IPeak>> spectrumFilterFunction;
    private final IComparisonPredicate<ICluster> clusterComparisonPredicate;

    private int currentMZAsInt;
    private INumberOfComparisonAssessor numberOfComparisonAssessor;

    public GreedyIncrementalClusteringEngine(ISimilarityChecker sck,
                                             Comparator<ICluster> scm,
                                             float windowSize,
                                             double clusteringPrecision,
                                             IFunction<List<IPeak>, List<IPeak>> spectrumFilterFunction,
                                             IComparisonPredicate<ICluster> clusterComparisonPredicate,
                                             INumberOfComparisonAssessor numberOfComparisonAssessor) {
        this.similarityChecker = sck;
        this.spectrumComparator = scm;
        this.windowSize = windowSize;
        // this change is performed so that a high threshold means a high clustering quality
        this.mixtureProbability = 1 - clusteringPrecision;
        this.spectrumFilterFunction = spectrumFilterFunction;
        this.clusterComparisonPredicate = clusterComparisonPredicate;
        this.numberOfComparisonAssessor = numberOfComparisonAssessor;

        try {
            this.cumulativeDistributionFunction = CumulativeDistributionFunctionFactory.getDefaultCumlativeDistributionFunctionForSimilarityMetric(sck.getClass());
        }
        catch(Exception e) {
            throw new IllegalStateException(e);
        }
    }

    public GreedyIncrementalClusteringEngine(ISimilarityChecker sck,
                                             Comparator<ICluster> scm,
                                             float windowSize,
                                             double clusteringPrecision,
                                             IFunction<List<IPeak>, List<IPeak>> spectrumFilterFunction,
                                             IComparisonPredicate<ICluster> clusterComparisonPredicate) {
        this(sck, scm, windowSize, clusteringPrecision, spectrumFilterFunction, clusterComparisonPredicate,
                Defaults.getNumberOfComparisonAssessor());
    }

    public GreedyIncrementalClusteringEngine(ISimilarityChecker sck,
                                             Comparator<ICluster> scm,
                                             float windowSize,
                                             double clusteringPrecision,
                                             IFunction<List<IPeak>, List<IPeak>> spectrumFilterFunction) {
        this(sck, scm, windowSize, clusteringPrecision, spectrumFilterFunction, null);
    }

    public double getWindowSize() {
        return windowSize;
    }


    public int getCurrentMZ() {
        return currentMZAsInt;
    }


    public void setCurrentMZ(final double pCurrentMZ) {
        int test = MZIntensityUtilities.mzToInt(pCurrentMZ);
        final int currentMZ = getCurrentMZ();
        if (currentMZ > test) {  // all ow roundoff error but not much
            double del = currentMZ - test;  // difference

            if (Math.abs(del) > MZIntensityUtilities.MZ_RESOLUTION * MZIntensityUtilities.SMALL_MZ_DIFFERENCE) {
                throw new IllegalStateException("mz values MUST be added in order - was "
                        + NumberUtilities.formatDouble(currentMZ, 3) + " new " +
                        NumberUtilities.formatDouble(pCurrentMZ, 3) + " del  " +
                        del
                );

            }

        }
        currentMZAsInt = test;
    }


    /**
     * Get clustered clusters
     */
    @Override
    public List<ICluster> getClusters() {
        final ArrayList<ICluster> ret = new ArrayList<>(clusters);
        Collections.sort(ret);
        return ret;
    }

    /**
     * add some clusters
     */
    @Override
    public void addClusters(ICluster... cluster) {
        // or use a WrappedIncrementalClusteringEngine
        throw new UnsupportedOperationException("Use addClusterIncremental instead or use a WrappedIncrementalClusteringEngine ");

    }

    /**
     * add one cluster and return any clusters which are too far in mz from further consideration
     *
     * @return !null Cluster
     */
    @Override
    public List<ICluster> addClusterIncremental(final ICluster added) {
        double precursorMz = added.getPrecursorMz();
        List<ICluster> clustersToremove = findClustersTooLow(precursorMz);
        // either add as an existing cluster if make a new cluster
        addToClusters(added);
        return clustersToremove;
    }

    /**
     * return a list of clusters whose mz is too low to merge with the current cluster
     * these are dropped and will be handled as never modifies this pass
     *
     * @param precursorMz new Mz
     * @return !null list of clusters to remove
     */
    protected List<ICluster> findClustersTooLow(double precursorMz) {
        setCurrentMZ(precursorMz); // also performs sanity check whether precursorMz is larger than currentMz

        double windowSize1 = getWindowSize();
        double lowestMZ = precursorMz - windowSize1;
        List<ICluster> clustersToremove = new ArrayList<>();
        List<ISpectrum> consensusSpectraToRemove = new ArrayList<>();

        for (int i = 0; i < clusters.size(); i++) {
            ICluster currentCluster = clusters.get(i);
            ISpectrum currentConsensusSpec = filteredConsensusSpectra.get(i);

            float testPrecursorMz = currentCluster.getPrecursorMz();
            if (lowestMZ > testPrecursorMz) {
                clustersToremove.add(currentCluster);
                consensusSpectraToRemove.add(currentConsensusSpec);
            }
        }
        if (!clustersToremove.isEmpty()) {
            clusters.removeAll(clustersToremove);
            filteredConsensusSpectra.removeAll(consensusSpectraToRemove);
        }

        return clustersToremove;
    }

    /**
     * this method is called by guaranteeClean to place any added clusters in play
     * for further clustering
     */
    protected void addToClusters(final ICluster clusterToAdd) {
        GreedySpectralCluster greedySpectralCluster = convertToGreedyCluster(clusterToAdd);

        // if there are no clusters yet, just save it
        if (clusters.isEmpty()) {
            clusters.add(greedySpectralCluster);
            filteredConsensusSpectra.add(filterSpectrum(greedySpectralCluster.getConsensusSpectrum()));
            return;
        }

        ISimilarityChecker sCheck = getSimilarityChecker();
        ISpectrum consensusSpectrumToAdd = clusterToAdd.getConsensusSpectrum();
        // always only compare the N highest peaks
        ISpectrum filteredConsensusSpectrumToAdd = filterSpectrum(consensusSpectrumToAdd);

        // add once an acceptable similarity score is found
        // this version does not look for the best match
        int nComparisons = numberOfComparisonAssessor.getNumberOfComparisons(clusterToAdd, clusters.size());

        for (int i = 0; i < clusters.size(); i++) {
            GreedySpectralCluster existingCluster = clusters.get(i);

            // apply the predicate if needed
            if (clusterComparisonPredicate != null) {
                if (!clusterComparisonPredicate.apply(clusterToAdd, existingCluster))
                    continue;
            }

            ISpectrum filteredConsensusSpectrum = filteredConsensusSpectra.get(i);

            double similarityScore = sCheck.assessSimilarity(filteredConsensusSpectrum, filteredConsensusSpectrumToAdd);

            if (cumulativeDistributionFunction.isSaveMatch(similarityScore, nComparisons, mixtureProbability)) {
                // use the originally passed cluster object for this, the greedy version is only used
                // to track comparison results and used if added internally

                // save the number of comparisons present when adding single spectra
                if (Defaults.isSaveDebugInformation()) {
                    if (clusterToAdd.getClusteredSpectraCount() == 1) {
                        clusterToAdd.getClusteredSpectra().get(0).setProperty(KnownProperties.MIN_COMPARISONS, String.valueOf(nComparisons));
                    }
                    if (existingCluster.getClusteredSpectraCount() == 1) {
                        existingCluster.getClusteredSpectra().get(0).setProperty(KnownProperties.MIN_COMPARISONS, String.valueOf(nComparisons));
                    }
                }

                // save the score if this is set
                if (Defaults.isSaveAddingScore()) {
                    // only save the score for single spectra
                    if (clusterToAdd.getClusteredSpectraCount() == 1) {
                        clusterToAdd.getClusteredSpectra().get(0).setProperty(
                                KnownProperties.ADDING_SCORE, String.valueOf(similarityScore));
                    }
                    if (existingCluster.getClusteredSpectraCount() == 1) {
                        existingCluster.getClusteredSpectra().get(0).setProperty(
                                KnownProperties.ADDING_SCORE, String.valueOf(similarityScore));
                    }
                }

                // preserve the id of the larger cluster
                if (clusterToAdd.getClusteredSpectraCount() > existingCluster.getClusteredSpectraCount())
                    existingCluster.setId(clusterToAdd.getId());

                // add to cluster
                existingCluster.addCluster(clusterToAdd);

                // update the existing consensus spectrum
                filteredConsensusSpectra.set(i, filterSpectrum(existingCluster.getConsensusSpectrum()));

                // since the cluster was added we're done
                return;
            }

            // save the comparison result for the next round of clustering
            greedySpectralCluster.saveComparisonResult(existingCluster.getId(), (float) similarityScore);
            existingCluster.saveComparisonResult(greedySpectralCluster.getId(), (float) similarityScore);
        }

        // since the cluster wasn't merged, add it as new
        clusters.add(greedySpectralCluster);
        // process the consensus spectrum
        ISpectrum filteredConsensusSpectrum = filterSpectrum(greedySpectralCluster.getConsensusSpectrum());
        filteredConsensusSpectra.add(filteredConsensusSpectrum);
    }

    private ISpectrum filterSpectrum(ISpectrum spectrumToFilter) {
        if (spectrumFilterFunction == null)
            return spectrumToFilter;

        ISpectrum filteredSpectrum;
        String nHighestPeaks = spectrumToFilter.getProperty(KnownProperties.N_HIGHEST_PEAKS);
        if (nHighestPeaks != null) {
            int highestPeaks = Integer.parseInt(nHighestPeaks);
            filteredSpectrum = spectrumToFilter.getHighestNPeaks(highestPeaks);
        }
        else {
            filteredSpectrum = new Spectrum(spectrumToFilter, spectrumFilterFunction.apply(spectrumToFilter.getPeaks()));
            filteredSpectrum.setProperty(KnownProperties.N_HIGHEST_PEAKS, String.valueOf(filteredSpectrum.getPeaks().size()));
        }

        return filteredSpectrum;
    }

    private GreedySpectralCluster convertToGreedyCluster(ICluster cluster) {
        // change the clusterToAdd to a 'GreedyCluster'
        GreedySpectralCluster greedyCluster;
        if (GreedySpectralCluster.class.isInstance(cluster)) {
            greedyCluster = (GreedySpectralCluster) cluster;
        }
        else {
            greedyCluster = new GreedySpectralCluster(cluster);
        }

        return greedyCluster;
    }

    /**
     * clusters are merged in the internal collection
     *
     * @return true is  anything happened
     */
    @Override
    public boolean processClusters() {
        throw new UnsupportedOperationException("Don\'t do this using an IncrementalClusteringEngine use a WrappedIncrementalClusteringEngine"); // ToDo
    }

    /**
     * used to expose internals for overriding classes only
     *
     * @return
     */
    @Override
    public ISimilarityChecker getSimilarityChecker() {
        return similarityChecker;
    }

    /**
     * Get similarity threshold used
     *
     * @return
     */
    @Override
    public double getSimilarityThreshold() {
        return 1 - mixtureProbability;
    }

    /**
     * used to expose internals for overridimg classes only
     *
     * @return
     */
    @SuppressWarnings("UnusedDeclaration")
    protected Comparator<ICluster> getSpectrumComparator() {
        return spectrumComparator;
    }

    /**
     * allow engines to be named
     *
     * @return
     */
    @Override
    public String toString() {
        int nClusters = size();
        final String name = this.getClass().getName();
        if (name != null)
            return name + " with " + nClusters;
        return super.toString();
    }

    /**
     * total number of clusters including queued clustersToAdd
     *
     * @return
     */
    @Override
    public int size() {
        return clusters.size();
    }
}

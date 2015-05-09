package uk.ac.ebi.pride.spectracluster.engine;

import uk.ac.ebi.pride.spectracluster.cdf.CumulativeDistributionFunction;
import uk.ac.ebi.pride.spectracluster.cdf.CumulativeDistributionFunctionFactory;
import uk.ac.ebi.pride.spectracluster.cluster.GreedySpectralCluster;
import uk.ac.ebi.pride.spectracluster.cluster.ICluster;
import uk.ac.ebi.pride.spectracluster.cluster.SpectralCluster;
import uk.ac.ebi.pride.spectracluster.consensus.GreedyConsensusSpectrum;
import uk.ac.ebi.pride.spectracluster.consensus.IConsensusSpectrumBuilder;
import uk.ac.ebi.pride.spectracluster.similarity.ISimilarityChecker;
import uk.ac.ebi.pride.spectracluster.spectrum.IPeak;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.spectrum.KnownProperties;
import uk.ac.ebi.pride.spectracluster.spectrum.Spectrum;
import uk.ac.ebi.pride.spectracluster.util.Defaults;
import uk.ac.ebi.pride.spectracluster.util.MZIntensityUtilities;
import uk.ac.ebi.pride.spectracluster.util.NumberUtilities;
import uk.ac.ebi.pride.spectracluster.util.function.IFunction;

import javax.annotation.Nonnull;
import java.util.*;

/**
 * uk.ac.ebi.pride.spectracluster.engine.IncrementalClusteringEngine
 * a version of a clustering engine in which spectra are added incrementally and
 * clusters are shed when they are too far to use
 * <p/>
 * <p/>
 * User: Steve
 * Date: 7/5/13
 */
public class GreedyIncrementalClusteringEngine implements IIncrementalClusteringEngine {
    private final List<GreedySpectralCluster> clusters = new ArrayList<GreedySpectralCluster>();
    private final ISimilarityChecker similarityChecker;
    private final Comparator<ICluster> spectrumComparator;
    private final double windowSize;
    private final double mixtureProbability;
    private final CumulativeDistributionFunction cumulativeDistributionFunction;
    private final IFunction<List<IPeak>, List<IPeak>> spectrumFilterFunction;

    private int currentMZAsInt;

    public GreedyIncrementalClusteringEngine(ISimilarityChecker sck,
                                             Comparator<ICluster> scm,
                                             float windowSize,
                                             double clusteringPrecision,
                                             IFunction<List<IPeak>, List<IPeak>> spectrumFilterFunction) {
        this.similarityChecker = sck;
        this.spectrumComparator = scm;
        this.windowSize = windowSize;
        // this change is performed so that a high threshold means a high clustering quality
        this.mixtureProbability = 1 - clusteringPrecision;
        this.spectrumFilterFunction = spectrumFilterFunction;
        try {
            this.cumulativeDistributionFunction = CumulativeDistributionFunctionFactory.getCumulativeDistributionFunctionForSimilarityMetric(sck.getClass());
        }
        catch(Exception e) {
            throw new IllegalStateException(e);
        }
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
        final ArrayList<ICluster> ret = new ArrayList<ICluster>(clusters);
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
        List<ICluster> clustersToremove = new ArrayList<ICluster>();
        for (ICluster test : clusters) {
            float testPrecursorMz = test.getPrecursorMz();
            if (lowestMZ > testPrecursorMz) {
                clustersToremove.add(test);
            }
        }
        if (!clustersToremove.isEmpty())
            clusters.removeAll(clustersToremove);

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
            return;
        }

        ISimilarityChecker sCheck = getSimilarityChecker();
        ISpectrum consensusSpectrumToAdd = clusterToAdd.getConsensusSpectrum();
        ISpectrum filteredConsensusSpectrumToAdd = filterSpectrum(consensusSpectrumToAdd);

        // add once an acceptable similarity score is found
        // this version does not look for the best match
        int nComparisons = 0;

        for (GreedySpectralCluster existingCluster : clusters) {
            ISpectrum consensusSpectrum = existingCluster.getConsensusSpectrum();
            ISpectrum filteredConsensusSpectrum = filterSpectrum(consensusSpectrum);

            double similarityScore = sCheck.assessSimilarity(filteredConsensusSpectrum, filteredConsensusSpectrumToAdd);
            nComparisons++;

            if (cumulativeDistributionFunction.isSaveMatch(similarityScore, nComparisons, mixtureProbability)) {
                // use the originally passed cluster object for this, the greedy version is only used
                // to track comparison results and used if added internally

                // preserve the id of the larger cluster
                if (clusterToAdd.getClusteredSpectraCount() > existingCluster.getClusteredSpectraCount())
                    existingCluster.setId(clusterToAdd.getId());

                // add to cluster
                existingCluster.addCluster(clusterToAdd);

                // since the cluster was added we're done
                return;
            }

            // save the comparison result for the next round of clustering
            greedySpectralCluster.saveComparisonResult(existingCluster.getId(), (float) similarityScore);
            existingCluster.saveComparisonResult(greedySpectralCluster.getId(), (float) similarityScore);
        }

        // since the cluster wasn't merged, add it as new
        clusters.add(greedySpectralCluster);
    }

    private ISpectrum filterSpectrum(ISpectrum spectrumToFilter) {
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

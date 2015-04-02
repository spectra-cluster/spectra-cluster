package uk.ac.ebi.pride.spectracluster.engine;

import uk.ac.ebi.pride.spectracluster.cluster.ICluster;
import uk.ac.ebi.pride.spectracluster.cluster.SpectralCluster;
import uk.ac.ebi.pride.spectracluster.similarity.ISimilarityChecker;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.util.Defaults;
import uk.ac.ebi.pride.spectracluster.util.MZIntensityUtilities;
import uk.ac.ebi.pride.spectracluster.util.NumberUtilities;

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
public class IncrementalClusteringEngine implements IIncrementalClusteringEngine {
    /**
     * Defines the proportion of spectra that must be shared between two clusters
     * to define these two clusters as identical. Comparison is performed based
     * on spectra ids.
     * Setting this value to 1 only merges clusters that contain the exact same
     * spectra.
     */
    public static final double PROPORTION_SHARED_SPECTRA_FOR_IDENTICAL = 1;

    /**
     * These are mainly for debugging
     */
    public static int numberNotMerge = 0;
    public static int numberGoodMerge = 0;
    public static int numberLessGoodMerge = 0;

    private final List<ICluster> clusters = new ArrayList<ICluster>();
    private final ISimilarityChecker similarityChecker;
    private final Comparator<ICluster> spectrumComparator;
    private final double windowSize;
    private final double similarityThreshold;
    private int currentMZAsInt;

    public IncrementalClusteringEngine(ISimilarityChecker sck,
                                       Comparator<ICluster> scm,
                                       float windowSize,
                                       double similarityThreshold) {
        this.similarityChecker = sck;
        this.spectrumComparator = scm;
        this.windowSize = windowSize;
        this.similarityThreshold = similarityThreshold;
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
    // TODO JG discuss whether parameter precursorMz is sensible or if the currentPrecurorMz should be used
    protected List<ICluster> findClustersTooLow(double precursorMz) {
        setCurrentMZ(precursorMz); // also performs sanity check whether precursorMz is larger than currentMz

        double windowSize1 = getWindowSize();
        double lowestMZ = precursorMz - windowSize1;
        List<ICluster> clustersToremove = new ArrayList<ICluster>();
        List<ICluster> myClusters = internalGetClusters();
        for (ICluster test : myClusters) {
            float testPrecursorMz = test.getPrecursorMz();
            if (lowestMZ > testPrecursorMz) {
                clustersToremove.add(test);
            }
        }
        if (!clustersToremove.isEmpty())
            internalGetClusters().removeAll(clustersToremove);   // might break hear

        return clustersToremove;

    }


    /**
     * this method is called by guaranteeClean to place any added clusters in play
     * for further clustering
     */
    protected void addToClusters(final ICluster clusterToAdd) {
        List<ICluster> myClusters = internalGetClusters();
        if (myClusters.isEmpty()) {   // no checks just add
            myClusters.add(new SpectralCluster(clusterToAdd, Defaults.getDefaultConsensusSpectrumBuilder()));
            numberNotMerge++;
            return;
        }

        //  handle the case of subsets and almost subsets
        if (handleFullContainment(clusterToAdd))
            return; // no need to add we are contained


        ISimilarityChecker sCheck = getSimilarityChecker();
        List<ISpectrum> clusteredSpectra1 = clusterToAdd.getClusteredSpectra();

        double highestSimilarityScore = 0;
        ICluster mostSimilarCluster = null;

        ISpectrum consensusSpectrum1 = clusterToAdd.getConsensusSpectrum();  // sub spectra are really only one spectrum clusters
        // find the cluster with the highest similarity score
        for (ICluster cluster : myClusters) {
            ISpectrum consensusSpectrum = cluster.getConsensusSpectrum();

            double similarityScore = sCheck.assessSimilarity(consensusSpectrum, consensusSpectrum1);

            if (similarityScore > highestSimilarityScore && similarityScore >= similarityThreshold) {
                highestSimilarityScore = similarityScore;
                mostSimilarCluster = cluster;
            }
        }

        if (mostSimilarCluster != null) {
            // add to cluster
            ISpectrum[] clusteredSpectra = new ISpectrum[clusteredSpectra1.size()];
            final ISpectrum[] merged = clusteredSpectra1.toArray(clusteredSpectra);
            mostSimilarCluster.addSpectra(merged);

            // Preserve the cluster id from the bigger cluster, in terms of number of spectra
            // This is used to facilitate incremental clustering
            if (mostSimilarCluster.getClusteredSpectraCount() < clusterToAdd.getClusteredSpectraCount()) {
                mostSimilarCluster.setId(clusterToAdd.getId());
            }

            numberGoodMerge++;
        } else {
            // create a new cluster
            myClusters.add(new SpectralCluster(clusterToAdd, Defaults.getDefaultConsensusSpectrumBuilder()));
            numberNotMerge++;
        }
    }

    /**
     * Merges clusters if they share more then PROPORTION_SHARED_SPECTRA_FOR_IDENTICAL spectra. The
     * proportion is calculated relative to the smaller cluster's size.
     *
     * @param clusterToAdd
     * @return true is we fully replace a cluster with a larger or find this fully contained
     */
    protected boolean handleFullContainment(final ICluster clusterToAdd) {
        final List<ICluster> myclusters = internalGetClusters();
        ICluster toReplace = null;
        double maxProportionSharedSpectra = Double.MIN_VALUE;
        for (ICluster myCluster : myclusters) {
            double proportionSharedSpectra = getProportionSharedSpectraIds(clusterToAdd, myCluster);

            if (proportionSharedSpectra > maxProportionSharedSpectra) {
                maxProportionSharedSpectra = proportionSharedSpectra;
                toReplace = myCluster;
                if (proportionSharedSpectra == 1)   // is full subset
                    break;
            }
        }

        if (maxProportionSharedSpectra >= PROPORTION_SHARED_SPECTRA_FOR_IDENTICAL) {
            mergeIntoCluster(clusterToAdd, toReplace);
            return true; // done
        }
        return false;
    }

    /**
     * return common spectra ids
     *
     * @param c2 cluster2
     * @return set of common ids
     */
    public static
    @Nonnull
    Set<String> getSharedSpectraIds(@Nonnull final Set<String> firstIds, @Nonnull final ICluster c2) {
        Set<String> ret = new HashSet<String>(firstIds);
        ret.retainAll(c2.getSpectralIds());
        return ret;
    }

    public static
    @Nonnull
    double getProportionSharedSpectraIds(@Nonnull final ICluster cluster1, @Nonnull final ICluster cluster2) {
        int sharedSpectraIds = getSharedSpectraIds(cluster1.getSpectralIds(), cluster2).size();

        int minSize = Math.min(cluster1.getClusteredSpectraCount(), cluster2.getClusteredSpectraCount());

        return (double) sharedSpectraIds / minSize;
    }

    protected void mergeIntoCluster(final ICluster mergeFrom, final ICluster mergeInto) {
        List<ISpectrum> clusteredSpectra1 = mergeFrom.getClusteredSpectra();
        ISpectrum[] clusteredSpectra = new ISpectrum[clusteredSpectra1.size()];
        final ISpectrum[] merged = clusteredSpectra1.toArray(clusteredSpectra);
        mergeInto.addSpectra(merged);
        numberLessGoodMerge++;
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
    protected List<ICluster> internalGetClusters() {
        return clusters;
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
        return similarityThreshold;
    }

    /**
     * used to expose internals for overriding classes only
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

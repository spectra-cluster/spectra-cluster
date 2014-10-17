package uk.ac.ebi.pride.spectracluster.engine;

import uk.ac.ebi.pride.spectracluster.cluster.ICluster;
import uk.ac.ebi.pride.spectracluster.cluster.SpectralCluster;
import uk.ac.ebi.pride.spectracluster.similarity.ISimilarityChecker;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.util.ClusterUtilities;
import uk.ac.ebi.pride.spectracluster.util.Defaults;
import uk.ac.ebi.pride.spectracluster.util.MZIntensityUtilities;
import uk.ac.ebi.pride.spectracluster.util.NumberUtilities;

import javax.annotation.Nonnull;
import java.util.*;

/**
 * uk.ac.ebi.pride.spectracluster.engine.IncrementalClusteringEngine
 * a version of a clustering enging in which spectra are added incrementatlly and
 * clusters are shed when they are too far to use
 * <p/>
 * <p/>
 * User: Steve
 * Date: 7/5/13
 */
public class IncrementalClusteringEngine implements IIncrementalClusteringEngine {


    public static final double MINIMUM_SIMILARITY_SCORE_FOR_OVERLAP = 0.2;
    public static final double BONUS_PER_OVERLAP = 0.05;
    // with 1 we merge all true subclusters but not others
    public static final double MINIMUM_MERGE_SCORE = 1; // 0.5;


    /**
     * These are mainly for debugging
     */
    public static int numberOverlap = 0;
    public static int numberNotMerge = 0;
    public static int numberGoodMerge = 0;
    public static int numberLessGoodMerge = 0;
    public static int numberReAsssigned = 0;

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
     * SLewis - I think a guarantee that they are sorted by MZ is useful
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
        double defaultThreshold1 = getWindowSize();
        double lowestMZ = precursorMz - defaultThreshold1;
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


        setCurrentMZ(precursorMz);
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
        if (handleFullContainment(clusterToAdd))    // TODO JG this should not be necessary. Clusters with subsets should be similar enough to be merged by the clustering step anyway
            return; // no need to add we are contained


        ISimilarityChecker sCheck = getSimilarityChecker();
        List<ISpectrum> clusteredSpectra1 = clusterToAdd.getClusteredSpectra();

        ICluster bestMatch = null;
        double highestSimilarityScore = 0;

        ICluster mostSimilarCluster = null;
        ISpectrum consensusSpectrum1 = clusterToAdd.getConsensusSpectrum();  // subspectra are really only one spectrum clusters
        // find the cluster with the highest similarity score
        for (ICluster cluster : myClusters) {
            ISpectrum consensusSpectrum = cluster.getConsensusSpectrum();

            double similarityScore = sCheck.assessSimilarity(consensusSpectrum, consensusSpectrum1);

            if (similarityScore > highestSimilarityScore) {
                highestSimilarityScore = similarityScore;
                bestMatch = cluster; //    track good but not great
                if (similarityScore >= similarityThreshold) {
                    mostSimilarCluster = bestMatch;
                }

            }
        }


        // add to cluster
        if (mostSimilarCluster != null) {
            ISpectrum[] clusteredSpectra = new ISpectrum[clusteredSpectra1.size()];
            final ISpectrum[] merged = clusteredSpectra1.toArray(clusteredSpectra);
            mostSimilarCluster.addSpectra(merged);
            numberGoodMerge++;
            return;
        }

        // maybe a lot of overlap here // TODO JG evaluate if this approach is sensible
        /*
        if (bestMatch != null) {
            if (handlePotentialOverlap(clusterToAdd, bestMatch, highestSimilarityScore))
                return;
        }
        */
        myClusters.add(new SpectralCluster(clusterToAdd, Defaults.getDefaultConsensusSpectrumBuilder()));
        numberNotMerge++;
    }

    /**
     * figure out is
     *
     * @param clusterToAdd
     * @return true is we fully replace a cluster with a larger or find this fully contained
     */
    protected boolean handleFullContainment(final ICluster clusterToAdd) {
        final List<ICluster> myclusters = internalGetClusters();
        ICluster toReplace = null;
        double bestSimilarity = Double.MIN_VALUE;
        for (ICluster myCluster : myclusters) {
            double score = ClusterUtilities.clusterFullyContainsScore(myCluster, clusterToAdd);
            if (score > bestSimilarity) {
                bestSimilarity = score;
                toReplace = myCluster;
                if (score == 1)   // is full subset
                    break;
            }
        }

        if (bestSimilarity >= MINIMUM_MERGE_SCORE) {
            mergeIntoCluster(clusterToAdd, toReplace);
            return true; // done
        }
        return false;
    }


    /**
     * if there are overlapping spectra among the current cluster and the best match
     * then  firure out what is best
     *
     * @param cluster1
     * @param cluster2
     * @return
     */
    @Deprecated // TODO JG Reevaluate the usage of this function
    protected boolean handlePotentialOverlap(final ICluster cluster1, final ICluster cluster2, double highestSimilarityScore) {
        if (highestSimilarityScore < MINIMUM_SIMILARITY_SCORE_FOR_OVERLAP)
            return false;     // we did nothing
        Set<String> ids = cluster1.getSpectralIds();
        Set<String> best = cluster2.getSpectralIds();
        Set<String> spectraOverlap = getSpectraOverlap(ids, cluster2);
        int numberOverlap = spectraOverlap.size();
        if (numberOverlap == 0)
            return false; // no overlap
        int minClusterSize = Math.min(best.size(), ids.size());

        // of a lot of overlap then force a merge
        if (numberOverlap >= minClusterSize / 2) {  // enough overlap then merge
            mergeIntoCluster(cluster1, cluster2);
            return true;
        }
        // allow a bonus for overlap
        ISimilarityChecker sCheck = getSimilarityChecker();
        if (highestSimilarityScore + BONUS_PER_OVERLAP > similarityThreshold) {
            mergeIntoCluster(cluster1, cluster2);
            return true;

        }


        // force overlappping spectra into the best cluster
        return assignOverlapsToBestCluster(cluster1, cluster2, spectraOverlap);


    }

    /**
     * return common spectra ids
     *
     * @param c2 cluster2
     * @return set of common ids
     */
    public static
    @Nonnull
    @Deprecated
    // TODO JG function highly similar to ClusterUtilities::clusterFullyContainsScore
    Set<String> getSpectraOverlap(@Nonnull final Set<String> firstIds, @Nonnull final ICluster c2) {
        Set<String> ret = new HashSet<String>(firstIds);
        ret.retainAll(c2.getSpectralIds());
        return ret;
    }

    protected void mergeIntoCluster(final ICluster mergeFrom, final ICluster mergeInto) {
        List<ISpectrum> clusteredSpectra1 = mergeFrom.getClusteredSpectra();
        ISpectrum[] clusteredSpectra = new ISpectrum[clusteredSpectra1.size()];
        final ISpectrum[] merged = clusteredSpectra1.toArray(clusteredSpectra);
        mergeInto.addSpectra(merged);
        numberLessGoodMerge++;
    }

    /**
     * this assigns
     *
     * @param cluster1
     * @param cluster2
     * @return
     */
    @Deprecated // TODO JG Clusters that overlap to a significant portion should be merged completely and not fragmented
    protected boolean assignOverlapsToBestCluster(final ICluster cluster1, final ICluster cluster2, Set<String> spectraOverlap) {
        List<ISpectrum> clusteredSpectra1 = cluster1.getClusteredSpectra();
        // I am not sure here but I think we let duplicates move to the best cluster
        ISimilarityChecker sCheck = getSimilarityChecker();
        ISpectrum cs1 = cluster1.getConsensusSpectrum();  // subspectra are really only one spectrum clusters
        ISpectrum cs2 = cluster2.getConsensusSpectrum();  // subspectra are really only one spectrum clusters
        for (ISpectrum spc : clusteredSpectra1) {
            if (!spectraOverlap.contains(spc.getId()))
                continue; // not an overlapping spectrum
            // choose the better cluster
            double ss1 = sCheck.assessSimilarity(cs1, spc);
            double ss2 = sCheck.assessSimilarity(cs2, spc);
            if (ss1 < ss2)
                cluster1.removeSpectra(spc);
            else
                cluster2.removeSpectra(spc);
            numberReAsssigned++;
        }

        return false;
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

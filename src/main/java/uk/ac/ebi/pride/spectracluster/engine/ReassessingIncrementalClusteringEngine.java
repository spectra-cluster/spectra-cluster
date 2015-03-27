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
 * ReassessingIncrementalClusteringEngine
 * This incremental clustering engine also clusters spectra incrementally and
 * sheds clusters when they are too far to use (based on precursor m/z). The main
 * difference is that it also checks all clusters for non-fitting spectra in case they
 * were changed.
 * @author Johannes Griss
  */
public class ReassessingIncrementalClusteringEngine implements IIncrementalClusteringEngine {
    /**
     * Defines the proportion of spectra that must be shared between two clusters
     * to define these two clusters as identical. Comparison is performed based
     * on spectra ids.
     * Setting this value to 1 only merges clusters that contain the exact same
     * spectra.
     */
    public static final double PROPORTION_SHARED_SPECTRA_FOR_IDENTICAL = 1;

    public static final String CLUSTER_PROPERTY_IS_DIRTY = "is_dirty";

    /**
     * These are mainly for debugging
     */
    public static int numberNotMerge = 0;
    public static int numberGoodMerge = 0;
    public static int numberLessGoodMerge = 0;

    public static final double DEFAULT_RETAIN_THRESHOLD = Defaults.getRetainThreshold();

    private final List<ICluster> clusters = new ArrayList<ICluster>();
    private final ISimilarityChecker similarityChecker;
    private final Comparator<ICluster> spectrumComparator;
    private final double windowSize;
    private final double similarityThreshold;
    private final double similarityRetainThreshold;
    private int currentMZAsInt;

    public ReassessingIncrementalClusteringEngine(ISimilarityChecker sck,
                                                  Comparator<ICluster> scm,
                                                  float windowSize,
                                                  double similarityThreshold) {
        this.similarityChecker = sck;
        this.spectrumComparator = scm;
        this.windowSize = windowSize;
        this.similarityThreshold = similarityThreshold;
        this.similarityRetainThreshold = DEFAULT_RETAIN_THRESHOLD;
    }

    public ReassessingIncrementalClusteringEngine(ISimilarityChecker similarityChecker, Comparator<ICluster> spectrumComparator, double windowSize, double similarityThreshold, double similarityRetainThreshold) {
        this.similarityChecker = similarityChecker;
        this.spectrumComparator = spectrumComparator;
        this.windowSize = windowSize;
        this.similarityThreshold = similarityThreshold;
        this.similarityRetainThreshold = similarityRetainThreshold;
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

        // remove any non-fitting spectra
        List<ICluster> nonFittingSpectra = removeNoneFittingSpectra();

        // try to add them again by merging all clusters
        internalGetClusters().addAll(nonFittingSpectra);
        mergeAllClusters();

        return clustersToremove;
    }

    /**
     * Merge the clusters only, all the merged cluster
     *
     * @return true if clusters have been merged
     */
    public boolean mergeAllClusters() {
        boolean modified = false;
        boolean toMerge = true;
        List<ICluster> myClusters = internalGetClusters();
        ISimilarityChecker sCheck = getSimilarityChecker();

        while (toMerge) {
            toMerge = false;
            List<ICluster> clustersToRemove = new ArrayList<ICluster>();
            for (int i = 0; i < myClusters.size(); i++) {
                for (int j = i + 1; j < myClusters.size(); j++) {
                    ICluster clusterI = myClusters.get(i);
                    ICluster clusterJ = myClusters.get(j);
                    double similarityScore = sCheck.assessSimilarity(clusterI.getConsensusSpectrum(), clusterJ.getConsensusSpectrum());
                    if (similarityScore >= similarityThreshold) {
                        toMerge = true;
                        modified = true;
                        ISpectrum[] clusteredSpectra = new ISpectrum[clusterI.getClusteredSpectra().size()];
                        clusterJ.addSpectra(clusterI.getClusteredSpectra().toArray(clusteredSpectra));
                        clusterJ.setProperty(CLUSTER_PROPERTY_IS_DIRTY, ""); // mark the cluster dirty
                        clustersToRemove.add(clusterI);
                        break;
                    }
                }
            }
            clusters.removeAll(clustersToRemove);
        }

        return modified;
    }

    /**
     * Remove none fitting spectra and the none fitting spectra
     * back in as new clusters
     *
     * @return The spectra that were removed from the clusters
     */
    protected List<ICluster> removeNoneFittingSpectra() {
        List<ICluster> emptyClusters = new ArrayList<ICluster>(); // holder for any empty clusters
        List<ICluster> myClusters = internalGetClusters();
        List<ICluster> removedSpectra = new ArrayList<ICluster>(); // holder for new clusters

        for (ICluster cluster : myClusters) {
            // ignore any cluster that hasn't been changed
            if (cluster.getProperty(CLUSTER_PROPERTY_IS_DIRTY) == null)
                continue;

            List<ICluster> noneFittingSpectra = ClusterUtilities.findNoneFittingSpectra(cluster, similarityChecker, similarityRetainThreshold);
            if (!noneFittingSpectra.isEmpty()) {
                List<ISpectrum> holder = new ArrayList<ISpectrum>();
                for (ICluster nonFittingSpectrum : noneFittingSpectra) {
                    holder.addAll(nonFittingSpectrum.getClusteredSpectra());
                    removedSpectra.add(nonFittingSpectrum);
                }

                ISpectrum[] spectraToRemove = holder.toArray(new ISpectrum[holder.size()]);
                cluster.removeSpectra(spectraToRemove);

                if (cluster.getClusteredSpectraCount() == 0) {
                    emptyClusters.add(cluster); // nothing left remember this cluster
                }
                else {
                    cluster.setProperty(CLUSTER_PROPERTY_IS_DIRTY, null); // show that all spectra match the consensus spectrum
                }
            }
        }
        if (!emptyClusters.isEmpty())    // any empty clusters
            clusters.removeAll(emptyClusters);   // drop them

        return removedSpectra;
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


        ISimilarityChecker theSimilarityChecker = getSimilarityChecker();
        List<ISpectrum> clusteredSpectraToAdd = clusterToAdd.getClusteredSpectra();

        double highestSimilarityScore = 0;
        ICluster mostSimilarCluster = null;

        ISpectrum consensusSpectrum1 = clusterToAdd.getConsensusSpectrum();  // subspectra are really only one spectrum clusters
        // find the cluster with the highest similarity score
        for (ICluster cluster : myClusters) {
            ISpectrum consensusSpectrum = cluster.getConsensusSpectrum();

            double similarityScore = theSimilarityChecker.assessSimilarity(consensusSpectrum, consensusSpectrum1);

            if (similarityScore > highestSimilarityScore && similarityScore >= similarityThreshold) {
                highestSimilarityScore = similarityScore;
                mostSimilarCluster = cluster;
            }
        }

        if (mostSimilarCluster != null) {
            // add to cluster
            ISpectrum[] clusteredSpectra = new ISpectrum[clusteredSpectraToAdd.size()];
            final ISpectrum[] merged = clusteredSpectraToAdd.toArray(clusteredSpectra);
            mostSimilarCluster.addSpectra(merged);
            mostSimilarCluster.setProperty(CLUSTER_PROPERTY_IS_DIRTY, "");
            numberGoodMerge++;
        }
        else {
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

package uk.ac.ebi.pride.spectracluster.engine;

import uk.ac.ebi.pride.spectracluster.cluster.GreedySpectralCluster;
import uk.ac.ebi.pride.spectracluster.cluster.ICluster;
import uk.ac.ebi.pride.spectracluster.similarity.ISimilarityChecker;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.util.Defaults;
import uk.ac.ebi.pride.spectracluster.util.MZIntensityUtilities;
import uk.ac.ebi.pride.spectracluster.util.NumberUtilities;

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
public class SimilarClusterMergingEngine implements IIncrementalClusteringEngine {
    protected final List<ICluster> clusters = new ArrayList<ICluster>();
    protected final Comparator<ICluster> spectrumComparator;
    protected final double windowSize;
    protected final double requiredSharedSpectra;

    protected int currentMZAsInt;

    public SimilarClusterMergingEngine(Comparator<ICluster> scm,
                                       float windowSize,
                                       double requiredSharedSpectra) {
        this.spectrumComparator = scm;
        this.windowSize = windowSize;
        this.requiredSharedSpectra = requiredSharedSpectra;
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
        Set<String> spectraIdsToAdd = clusterToAdd.getSpectralIds();

        for (ICluster existingCluster : clusters) {
            Set<String> existingSpectraIds = existingCluster.getSpectralIds();
            double sharedSpectra = calculateSharedSpectra(spectraIdsToAdd, existingSpectraIds);

            if (sharedSpectra >= requiredSharedSpectra) {
                // ignore perfect duplicates
                if (sharedSpectra == 1)
                    return;

                if (clusterToAdd.storesPeakLists()) {
                    ISpectrum[] buffer = new ISpectrum[clusterToAdd.getClusteredSpectra().size()];
                    existingCluster.addSpectra(clusterToAdd.getClusteredSpectra().toArray(buffer));
                }
                else {
                    // it must be a greedy cluster
                    // TODO: this kind of merging may lead to incorrect results
                    if (!GreedySpectralCluster.class.isInstance(existingCluster)) {
                        throw new IllegalStateException("Cannot add greedy cluster to non-greedy cluster");
                    }

                    GreedySpectralCluster greedySpectralCluster = (GreedySpectralCluster) existingCluster;
                    greedySpectralCluster.addCluster(clusterToAdd);
                }
            }
        }

        // since the cluster wasn't merged, add it as new
        clusters.add(clusterToAdd);
    }

    protected double calculateSharedSpectra(Set<String> spectraIdsToAdd, Set<String> existingSpectraIds) {
        Set<String> sharedIds = new HashSet<String>(spectraIdsToAdd);
        sharedIds.retainAll(existingSpectraIds);

        int minSize = Math.min(spectraIdsToAdd.size(), existingSpectraIds.size());

        return (double) sharedIds.size() / minSize;
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
        return Defaults.getDefaultSimilarityChecker();
    }

    /**
     * Get similarity threshold used
     *
     * @return
     */
    @Override
    public double getSimilarityThreshold() {
        return requiredSharedSpectra;
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

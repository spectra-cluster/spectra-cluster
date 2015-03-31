package uk.ac.ebi.pride.spectracluster.engine;

import uk.ac.ebi.pride.spectracluster.cluster.ICluster;
import uk.ac.ebi.pride.spectracluster.cluster.SpectralCluster;
import uk.ac.ebi.pride.spectracluster.similarity.ISimilarityChecker;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.util.ClusterUtilities;
import uk.ac.ebi.pride.spectracluster.util.Defaults;

import java.util.*;

/**
 * Default implementation of the clustering engine
 *
 * @author Steve Lewis
 * @author Rui Wang
 * @version $Id$
 */
public class ClusteringEngine implements IClusteringEngine {

    private boolean dirty;
    private final List<ICluster> clusters = new ArrayList<ICluster>();
    private final List<ICluster> clustersToAdd = new ArrayList<ICluster>();
    private final ISimilarityChecker similarityChecker;
    private final Comparator<ICluster> spectrumComparator;
    /**
     * Similarity threshold above which spectra are added to a cluster
     */
    private final double similarityThreshold;
    /**
     * Similarity threshold below which spectra are removed from a cluster
     */
    private final double retainThreshold;

    public ClusteringEngine(ISimilarityChecker similarityChecker,
                            Comparator<ICluster> spectrumComparator,
                            double similarityThreshold,
                            double retainThreshold) {
        this.similarityChecker = similarityChecker;
        this.spectrumComparator = spectrumComparator;
        this.similarityThreshold = similarityThreshold;
        this.retainThreshold = retainThreshold;
    }

    protected void guaranteeClean() {
        if (isDirty()) {
            filterClustersToAdd();
            List<ICluster> myClustersToAdd = getClustersToAdd();
            Collections.sort(myClustersToAdd, getSpectrumComparator());
            addToClusters();
            setDirty(false);
        }
    }

    /**
     * Remove clusters which are size zero
     */
    protected void filterClustersToAdd() {
        List<ICluster> myClustersToAdd = getClustersToAdd();
        List<ICluster> l2 = new ArrayList<ICluster>();
        for (ICluster sc : myClustersToAdd) {
            if (sc.getClusteredSpectraCount() > 0)
                l2.add(sc);
        }
        myClustersToAdd.clear();
        myClustersToAdd.addAll(l2);
    }

    /**
     * Get clustered clusters
     * SLewis - I think a guarantee that they are sorted by MZ is useful
     */
    @Override
    public List<ICluster> getClusters() {
        guaranteeClean();
        final ArrayList<ICluster> ret = new ArrayList<ICluster>(clusters);
        Collections.sort(ret);
        return ret;
    }

    protected boolean isDirty() {
        return dirty;
    }

    protected void setDirty(boolean dirty) {
        this.dirty = dirty;
    }

    /**
     * add some clusters
     */
    @Override
    public void addClusters(ICluster... cluster) {
        List<ICluster> myClustersToAdd = getClustersToAdd();
        if (cluster != null) {
            myClustersToAdd.addAll(Arrays.asList(cluster));
            setDirty(true);
        }

    }

    /**
     * this method is called by guaranteeClean to place any added clusters in play
     * for further clustering
     */
    protected void addToClusters() {
        List<ICluster> myClustersToAdd = getClustersToAdd();
        List<ICluster> myClusters = internalGetClusters();
        ISimilarityChecker sCheck = getSimilarityChecker();

        for (ICluster clusterToAdd : myClustersToAdd) {

            ICluster mostSimilarCluster = null;
            double highestSimilarityScore = 0;

            // find the cluster with the highest similarity score
            for (ICluster cluster : myClusters) {
                ISpectrum consensusSpectrum = cluster.getConsensusSpectrum();
                ISpectrum consensusSpectrum1 = clusterToAdd.getConsensusSpectrum();  // subspectra are really only one spectrum clusters

                double similarityScore = sCheck.assessSimilarity(consensusSpectrum, consensusSpectrum1);

                if (similarityScore >= similarityThreshold && similarityScore > highestSimilarityScore) {
                    highestSimilarityScore = similarityScore;
                    mostSimilarCluster = cluster;
                }
            }

            // add to cluster
            if (mostSimilarCluster != null) {
                ISpectrum[] clusteredSpectra = new ISpectrum[clusterToAdd.getClusteredSpectra().size()];
                mostSimilarCluster.addSpectra(clusterToAdd.getClusteredSpectra().toArray(clusteredSpectra));

                // Preserve the cluster id from the bigger cluster, in terms of number of spectra
                // This is used to facilitate incremental clustering
                if (mostSimilarCluster.getClusteredSpectraCount() < clusterToAdd.getClusteredSpectraCount()) {
                    mostSimilarCluster.setId(clusterToAdd.getId());
                }
            } else {
                myClusters.add(new SpectralCluster(clusterToAdd, Defaults.getDefaultConsensusSpectrumBuilder()));
            }
        }

        myClustersToAdd.clear();
    }

    /**
     * clusters are merged in the internal collection
     *
     * @return true is  anything happened
     */
    @Override
    public boolean processClusters() {
        guaranteeClean();

        if (size() < 2)
            return false; // nothing to do

        // merge clusters
        boolean merged = mergeAllClusters();

        // remove none fitting spectra
        boolean noneFittingSpectraFound = demergeNoneFittingSpectra();

        // set dirty state
        setDirty(noneFittingSpectraFound);

        return merged || noneFittingSpectraFound;
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
     * @return true if spectra have been removed
     */
    protected boolean demergeNoneFittingSpectra() {
        boolean noneFittingSpectraFound = false;

        List<ICluster> emptyClusters = new ArrayList<ICluster>(); // holder for any empty clusters
        List<ICluster> myClusters = internalGetClusters();

        for (ICluster cluster : myClusters) {
            List<ICluster> noneFittingSpectra = ClusterUtilities.findNoneFittingSpectra(cluster, similarityChecker, retainThreshold);
            if (!noneFittingSpectra.isEmpty()) {
                noneFittingSpectraFound = true;

                List<ISpectrum> holder = new ArrayList<ISpectrum>();
                for (ICluster removedCluster : noneFittingSpectra) {
                    holder.addAll(removedCluster.getClusteredSpectra());
                    clustersToAdd.add(removedCluster);
                }

                ISpectrum[] spectraToRemove = holder.toArray(new ISpectrum[holder.size()]);
                cluster.removeSpectra(spectraToRemove);

                if (cluster.getClusteredSpectraCount() == 0) {
                    emptyClusters.add(cluster); // nothing left remember this cluster
                }
            }
        }
        if (!emptyClusters.isEmpty())    // any empty clusters
            clusters.removeAll(emptyClusters);   // drop them

        return noneFittingSpectraFound;
    }

    /**
     * used to expose internals for overridig classes only
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
     * used to expose internals for overridig classes only
     *
     * @return
     */
    protected List<ICluster> getClustersToAdd() {
        return clustersToAdd;
    }

    /**
     * used to expose internals for overridig classes only
     *
     * @return
     */
    protected List<ICluster> internalGetClusters() {
        return clusters;
    }

    /**
     * used to expose internals for overridig classes only
     *
     * @return
     */
    protected Comparator<ICluster> getSpectrumComparator() {
        return spectrumComparator;
    }

    /**
     * total number of clusters including queued clustersToAdd
     *
     * @return
     */
    @Override
    public int size() {
        return clustersToAdd.size() + clusters.size();
    }
}

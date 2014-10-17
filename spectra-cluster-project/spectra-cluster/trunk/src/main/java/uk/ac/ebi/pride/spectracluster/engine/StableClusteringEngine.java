package uk.ac.ebi.pride.spectracluster.engine;

import uk.ac.ebi.pride.spectracluster.cluster.ICluster;
import uk.ac.ebi.pride.spectracluster.similarity.ISimilarityChecker;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;

import java.util.*;

/**
 * @author Steve Lewis
 * @author Rui Wang
 * @version $Id$
 *          <p/>
 *          todo: development
 */
@Deprecated
public class StableClusteringEngine implements IStableClusteringEngine {

    private final SortedSet<ICluster> unstableClusters = new TreeSet<ICluster>();

    private final ISimilarityChecker similarityChecker;

    private final double similarityThreshold;

    private boolean stableClusterProcessed;

    public StableClusteringEngine(ISimilarityChecker similarityChecker, double similarityThreshold) {
        this.similarityChecker = similarityChecker;
        this.similarityThreshold = similarityThreshold;
    }

    /**
     * add one cluster and return any clusters which are too far in mz from further consideration
     * NOTE clusters MUST be added in ascending MZ order
     *
     * @param added !null cluster to add
     * @return !null list of clusters not far enough away they will no longer change
     */
    @Override
    public Collection<ICluster> addClusterIncremental(final ICluster added) {
        throw new UnsupportedOperationException("This Should NEVER be Called");
    }

    /**
     * add some clusters
     *
     * @param cluster
     */
    @Override
    public void addClusters(final ICluster... cluster) {
        throw new UnsupportedOperationException("This Should NEVER be Called");
    }

    /**
     * clusters are merged in the internal collection
     *
     * @return true is  anything happened
     */
    @Override
    public boolean processClusters() {
        throw new UnsupportedOperationException("This Should NEVER be Called");
    }

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
     * total number of clusters including queued clustersToAdd
     *
     * @return
     */
    @Override
    public int size() {
        return unstableClusters.size();
    }

    @Override
    public void addUnstableCluster(ICluster unstableCluster) {
        if (stableClusterProcessed) {
            // let this slide
            // throw new IllegalStateException("Adding unstable cluster after processing stable clusters");
        }

        unstableClusters.add(unstableCluster);
    }

    @Override
    public void processStableCluster(ICluster stableCluster) {
        stableClusterProcessed = true;

        Set<ICluster> emptyUnstableClusters = new HashSet<ICluster>();
        for (ICluster unstableCluster : unstableClusters) {
            boolean empty = mergeUnstableCluster(stableCluster, unstableCluster);
            if (empty) {
                emptyUnstableClusters.add(unstableCluster);
            }
        }

        unstableClusters.removeAll(emptyUnstableClusters);
    }

    private boolean mergeUnstableCluster(ICluster stableCluster, ICluster unstableCluster) {
        ISpectrum consensusSpectrum = stableCluster.getConsensusSpectrum();

        // find all the unstable spectra which can be merged into the stable cluster
        Set<ISpectrum> spectraToRemove = new HashSet<ISpectrum>();
        for (ISpectrum unstableSpectrum : unstableCluster.getClusteredSpectra()) {
            double similarity = similarityChecker.assessSimilarity(unstableSpectrum, consensusSpectrum);
            if (similarity >= similarityThreshold) {
                spectraToRemove.add(unstableSpectrum);
            }
        }

        // remove from unstable cluster
        ISpectrum[] arrayOfSpectraToRemove = spectraToRemove.toArray(new ISpectrum[spectraToRemove.size()]);
        unstableCluster.removeSpectra(arrayOfSpectraToRemove);

        // add into stable cluster
        stableCluster.addSpectra(arrayOfSpectraToRemove);

        return unstableCluster.getClusteredSpectraCount() == 0;
    }

    @Override
    public Collection<ICluster> getClusters() {
        return new ArrayList<ICluster>(unstableClusters);
    }

    public int getNumberOfUnstableSpectra() {
        int count = 0;

        for (ICluster unstableCluster : unstableClusters) {
            count += unstableCluster.getClusteredSpectraCount();
        }

        return count;
    }
}

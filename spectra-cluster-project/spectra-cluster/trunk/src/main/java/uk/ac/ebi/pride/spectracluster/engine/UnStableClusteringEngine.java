package uk.ac.ebi.pride.spectracluster.engine;

import uk.ac.ebi.pride.spectracluster.cluster.ICluster;
import uk.ac.ebi.pride.spectracluster.similarity.ISimilarityChecker;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.Set;

/**
 * @author Steve Lewis
 * @author Rui Wang
 * @version $Id$
 *          <p/>
 *          todo: development
 */
@Deprecated
public class UnStableClusteringEngine implements IUnStableClusteringEngine {

    private final Set<ICluster> stableClusters = new HashSet<ICluster>();

    private final ISimilarityChecker similarityChecker;

    private boolean unStableClusterProcessed;

    private double similarityThreshold;


    public UnStableClusteringEngine(ISimilarityChecker similarityChecker, double similarityThreshold) {
        this.similarityChecker = similarityChecker;
        this.similarityThreshold = similarityThreshold;
    }

    @Override
    public void addStableCluster(ICluster stableCluster) {
        if (unStableClusterProcessed) {
            throw new IllegalStateException("Adding stable cluster after processing unstable clusters");
        }

        stableClusters.add(stableCluster);
    }


    /**
     * try to move spectra to the stable cluster
     *
     * @return true if changed
     */
    @Override
    public boolean processUnStableCluster(ICluster unStableCluster) {
        unStableClusterProcessed = true;

        int startSpectraCount = unStableCluster.getClusteredSpectraCount();
        for (ICluster stableCluster : stableClusters) {
            boolean empty = mergeUnstableCluster(stableCluster, unStableCluster);
            if (empty)
                break;
        }
        return startSpectraCount != unStableCluster.getClusteredSpectraCount();

    }

    private boolean mergeUnstableCluster(ICluster stableCluster, ICluster unstableCluster) {
        ISpectrum consensusSpectrum = stableCluster.getConsensusSpectrum();

        // find all the unstable spectra which can be merged into the stable cluster
        HashSet<ISpectrum> spectraToRemove = new HashSet<ISpectrum>();
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

    public Collection<ICluster> getClusters() {
        return new ArrayList<ICluster>(stableClusters);
    }

    public int getNumberOfUnstableSpectra() {
        int count = 0;

        for (ICluster unstableCluster : stableClusters) {
            count += unstableCluster.getClusteredSpectraCount();
        }

        return count;
    }
}

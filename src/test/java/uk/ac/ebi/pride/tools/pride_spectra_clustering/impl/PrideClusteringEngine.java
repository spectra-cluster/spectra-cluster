package uk.ac.ebi.pride.tools.pride_spectra_clustering.impl;

import uk.ac.ebi.pride.spectracluster.cluster.*;
import uk.ac.ebi.pride.spectracluster.engine.*;
import uk.ac.ebi.pride.spectracluster.similarity.*;
import uk.ac.ebi.pride.tools.pride_spectra_clustering.*;
import uk.ac.ebi.pride.tools.pride_spectra_clustering.util.*;

import java.util.*;

/**
 * Implement a clustering Engine using the original johannes code
 * NOT Implemented yet
 *
 * @author Steve Lewis
 * @version $Id$
 */
public class PrideClusteringEngine implements IClusteringEngine {

    private final List<ICluster> clusters = new ArrayList<>();
    private List<SpectraCluster> clustersFound = null;
    private final List<ClusteringSpectrum> addedSpectra = new ArrayList<>();
    private final SpectraClustering clustering = new FrankEtAlClustering();

    public PrideClusteringEngine() {
        clustering.setClusteringRounds(2);
        clustering.setSimilarityThreshold(0.7);
    }


    /**
     * Get clustered clusters
     * SLewis - I think a guarantee that they are sorted by MZ is useful
     */
    @Override
    public List<ICluster> getClusters() {
        //      guaranteeClean();
        return new ArrayList<>(clusters);
    }


    /**
     * add some clusters
     */
    @Override
    public void addClusters(ICluster... cluster) {
        if (cluster != null) {
            for (ICluster sc : cluster) {
                final SpectraCluster spectraCluster = Adapters.fromSpectraCluster(sc);
                final List<ClusteringSpectrum> spectra = spectraCluster.getSpectra();
                addedSpectra.addAll(spectra);
            }

        }
        clustersFound = null;

    }


    /**
     * clusters are merged in the internal collection
     *
     * @return true is  anything happened
     */
    @Override
    public boolean processClusters() {
        if (clustersFound != null)
            return false; // already done

        clusters.clear();
        clustersFound = clustering.clusterConvertedSpectra(addedSpectra);
        for (SpectraCluster cluster : clustersFound) {
            final ICluster spectralCluster = Adapters.fromSpectraCluster(cluster);
            clusters.add(spectralCluster);
        }

        return false; // we are done after one pass
    }

    @Override
    public ISimilarityChecker getSimilarityChecker() {
        throw new UnsupportedOperationException("Method call not supported in this class");
    }

    @Override
    public double getSimilarityThreshold() {
        return clustering.getSimilarityThreshold();
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

package uk.ac.ebi.pride.tools.pride_spectra_clustering;

import uk.ac.ebi.pride.tools.jmzreader.model.Spectrum;
import uk.ac.ebi.pride.tools.pride_spectra_clustering.util.ClusteringSpectrum;
import uk.ac.ebi.pride.tools.pride_spectra_clustering.util.SpectraCluster;

import java.util.List;

public interface SpectraClustering {
    /**
     * Clusters the passed spectra and returns
     * the resulting list of SpectraClusterS.
     *
     * @param spectra
     * @return
     */
    List<SpectraCluster> clusterSpectra(List<Spectrum> spectra);

    /**
     * Clusters the passed spectra and returns
     * the resulting list of SpectraClusterS.
     *
     * @param spectra
     * @return
     */
    List<SpectraCluster> clusterConvertedSpectra(List<ClusteringSpectrum> spectra);

    /**
     * Clusters the passed spectra base on already existing
     * cluster (ie. from a previous process) and returns
     * the resulting list of SpectraClusterS.
     *
     * @param spectra
     * @return
     */
    List<SpectraCluster> clusterConvertedSpectra(List<ClusteringSpectrum> spectra, List<SpectraCluster> existingCluster);

    void setSimilarityThreshold(Double threshold);

    Double getSimilarityThreshold();

    /**
     * Returns the number of times the clustering process
     * is repeated
     *
     * @return
     */
    int getClusteringRounds();

    /**
     * Sets the number of times the clustering process
     * is repeated
     *
     * @param rounds
     */
    void setClusteringRounds(int rounds);

    /**
     * Returns a (human readable) description of the clustering
     * process as well as the set parameters.
     *
     * @return
     */
    String getDescription();
}

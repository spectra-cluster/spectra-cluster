package uk.ac.ebi.pride.spectracluster.cluster;

import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;

/**
 * uk.ac.ebi.pride.spectracluster.cluster.ISpectrumHolder
 * generalize the concept holding spectra - ISpectralCluster can do
 * this but also new concensusSpectrumBuilder
 * User: Steve
 * Date: 7/10/13
 */
public interface ISpectrumHolder {

    /**
     * Add a array of spectrum to cluster
     */
    void addSpectra(ISpectrum... merged);


    /**
     * stable clusters do not support remove others do
     *
     * @return as above
     */
    boolean isRemoveSupported();

    /**
     * Remove an array of spectrum from cluster
     */
    void removeSpectra(ISpectrum... removed);

    /**
     * add a change listener
     * final to make sure this is not duplicated at multiple levels
     *
     * @param added non-null change listener
     */
    void addSpectrumHolderListener(SpectrumHolderListener added);

    /**
     * remove a change listener
     *
     * @param removed non-null change listener
     */
    void removeSpectrumHolderListener(SpectrumHolderListener removed);
}

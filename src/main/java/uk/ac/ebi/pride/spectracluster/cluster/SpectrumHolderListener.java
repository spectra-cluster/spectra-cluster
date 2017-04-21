package uk.ac.ebi.pride.spectracluster.cluster;

import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;

import java.io.*;

/**
 * uk.ac.ebi.pride.spectracluster.cluster.SpectrumHolderListener
 * User: Steve
 * Date: 7/10/13
 */
public interface SpectrumHolderListener extends Serializable {

    /**
     * handle notification of adding spectra
     *
     * @param holder !null holder
     * @param added  added spectra
     */
    void onSpectraAdd(ISpectrumHolder holder, ISpectrum... added);

    /**
     * handle notification of removing spectra
     *
     * @param holder  !null holder
     * @param removed removed spectra
     */
    void onSpectraRemove(ISpectrumHolder holder, ISpectrum... removed);

}

package uk.ac.ebi.pride.spectracluster.util;

import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;

import java.io.Serializable;

/**
 * uk.ac.ebi.pride.spectracluster.util.SpectrumCreateListener
 * make a listener for Spectrum reads
 * User: Steve
 * Date: 9/23/13
 */
public interface SpectrumCreateListener extends Serializable{

    /**
     * initialize reading - if reading happens once - say  from
     * one file all this may happen in the constructor
     */
    void onSpectrumStarted();

    /**
     * do something when a Spectrum is created or read
     *
     * @param Spectrum
     */
    void onSpectrumCreate(ISpectrum Spectrum);

    /**
     * do something when a Spectrum when the last Spectrum is read -
     * this may be after a file read is finished
     */
    void onSpectrumCreateFinished();


}

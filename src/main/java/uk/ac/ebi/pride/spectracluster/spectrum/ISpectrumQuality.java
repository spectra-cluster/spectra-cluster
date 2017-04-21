package uk.ac.ebi.pride.spectracluster.spectrum;

import uk.ac.ebi.pride.spectracluster.quality.*;

import java.io.*;

/**
 * Spectrum quality interface
 *
 * @author Rui Wang
 * @version $Id$
 */
public interface ISpectrumQuality extends Serializable {

    /**
     * Get the quality measure of a spectrum
     */
    double getQualityScore();


    /**
     * Get the quality scorer used
     *
     * @return quality scorer
     */
    IQualityScorer getQualityScorer();
}

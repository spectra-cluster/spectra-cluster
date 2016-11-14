package uk.ac.ebi.pride.spectracluster.quality;

import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.util.IAlgorithm;

/**
 * Calculate the quality score of a given spectrum and
 * return the quality as a double
 *
 * @author Rui Wang
 * @version $Id$
 */
public interface IQualityScorer extends IAlgorithm {
    /**
     * return the quality of the spectrum
     *
     * @param spectrum !null spectrum
     * @return quality &gt;= 0
     */
    double calculateQualityScore(ISpectrum spectrum);
}

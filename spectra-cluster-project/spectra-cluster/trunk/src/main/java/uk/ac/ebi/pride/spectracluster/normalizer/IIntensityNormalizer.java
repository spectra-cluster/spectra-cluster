package uk.ac.ebi.pride.spectracluster.normalizer;

import uk.ac.ebi.pride.spectracluster.spectrum.IPeak;
import uk.ac.ebi.pride.spectracluster.util.IAlgorithm;

import java.util.List;

/**
 * Normalizes a spectrum's intensities.
 *
 * @author Rui Wang
 */
public interface IIntensityNormalizer extends IAlgorithm {

    /**
     * return the value normalized to - especial;ly useful for total intensity normalization where
     * we may not weed to normalize
     *
     * @return as above
     */
    double getNormalizedValue();

    /**
     * normalize alist of peaks - all the dirty work is here
     *
     * @param peaks
     * @return
     */
    List<IPeak> normalizePeaks(List<IPeak> peaks);
}

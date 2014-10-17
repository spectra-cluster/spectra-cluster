package uk.ac.ebi.pride.spectracluster.filter;

import uk.ac.ebi.pride.spectracluster.spectrum.IPeak;

import java.util.List;

/**
 * IPeakFilter is an interface to filter a given list of peaks
 *
 * @author Rui Wang
 * @version $Id$
 */
public interface IPeakFilter {
    /**
     * a version of the filter that does nothing
     */
    public static final IPeakFilter NULL_FILTER = new IPeakFilter() {
        @Override public List<IPeak> filter(List<IPeak> peaks) {
                return peaks;
        }
    };
    /**
     * Filter a given list of peaks
     *
     * @param peaks given list of peaks
     * @return a list of filtered peaks
     */
    public List<IPeak> filter(List<IPeak> peaks);
}

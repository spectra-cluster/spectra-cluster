package uk.ac.ebi.pride.spectracluster.filter;

import uk.ac.ebi.pride.spectracluster.spectrum.IPeak;
import uk.ac.ebi.pride.spectracluster.util.comparator.PeakIntensityComparator;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

/**
 * uk.ac.ebi.pride.spectracluster.filter.HighestNPeakFilter
 * return the highest peaks up to  maxPeaks
 * default constructor uses 100
 *
 * @author Steve Lewis
 * @date 27/05/2014
 */
public class HighestNPeakFilter implements IPeakFilter {

    public static final int DEFAULT_MAX_PEAKS = 100;
    public static final Comparator<IPeak> INTENSITY_COMPARATOR = PeakIntensityComparator.INSTANCE;

    private final int maxPeaks;

    public HighestNPeakFilter(int maxPeaks) {
        this.maxPeaks = maxPeaks;
    }

    public HighestNPeakFilter() {
        this(DEFAULT_MAX_PEAKS);
    }


    /**
     * Filter a given list of peaks
     *
     * @param peaks given list of peaks
     * @return a list of filtered peaks
     */
    @Override
    public List<IPeak> filter(List<IPeak> peaks) {
        List<IPeak> byIntensity = new ArrayList<IPeak>(peaks);
        Collections.sort(byIntensity, INTENSITY_COMPARATOR);
        List<IPeak> ret = new ArrayList<IPeak>();
        for (IPeak iPeak : byIntensity) {
            ret.add(iPeak);
            if (ret.size() >= maxPeaks)
                break;
        }
        return ret;
    }
}

package uk.ac.ebi.pride.spectracluster.util.function.peak;

import uk.ac.ebi.pride.spectracluster.spectrum.IPeak;
import uk.ac.ebi.pride.spectracluster.spectrum.Peak;
import uk.ac.ebi.pride.spectracluster.util.comparator.PeakIntensityComparator;
import uk.ac.ebi.pride.spectracluster.util.function.IFunction;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

/**
 * Return the highest peaks up to  maxPeaks
 * default constructor uses 100
 *
 * @author Steve Lewis
 * @author Rui Wang
 * @version $Id$
 */
public class HighestNPeakFunction implements IFunction<List<IPeak>, List<IPeak>> {

    public static final int DEFAULT_MAX_PEAKS = 100;
    public static final Comparator<IPeak> INTENSITY_COMPARATOR = PeakIntensityComparator.INSTANCE;

    private final int maxPeaks;

    public HighestNPeakFunction() {
        this(DEFAULT_MAX_PEAKS);
    }

    public HighestNPeakFunction(int maxPeaks) {
        this.maxPeaks = maxPeaks;
    }

    @Override
    public List<IPeak> apply(List<IPeak> originalPeaks) {
        List<IPeak> byIntensity = new ArrayList<>(originalPeaks);
        byIntensity.sort(INTENSITY_COMPARATOR);
        List<IPeak> ret = new ArrayList<>();
        for (IPeak originalPeak : byIntensity) {
            ret.add(new Peak(originalPeak));
            if (ret.size() >= maxPeaks)
                break;
        }
        return ret;
    }
}

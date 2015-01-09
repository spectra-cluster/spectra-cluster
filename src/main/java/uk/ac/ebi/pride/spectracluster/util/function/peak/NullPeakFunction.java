package uk.ac.ebi.pride.spectracluster.util.function.peak;

import uk.ac.ebi.pride.spectracluster.spectrum.IPeak;
import uk.ac.ebi.pride.spectracluster.spectrum.Peak;
import uk.ac.ebi.pride.spectracluster.util.function.IFunction;

import java.util.ArrayList;
import java.util.List;

/**
 * Simply copy the original peaks
 *
 * @author Rui Wang
 * @version $Id$
 */
public class NullPeakFunction implements IFunction<List<IPeak>, List<IPeak>> {

    @Override
    public List<IPeak> apply(List<IPeak> originalPeaks) {
        List<IPeak> peaks = new ArrayList<IPeak>();
        for (IPeak originalPeak : originalPeaks) {
            peaks.add(new Peak(originalPeak));
        }

        return peaks;
    }
}

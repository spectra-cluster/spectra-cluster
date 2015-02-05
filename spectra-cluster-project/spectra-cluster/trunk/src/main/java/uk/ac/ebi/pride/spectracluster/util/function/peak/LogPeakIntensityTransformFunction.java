package uk.ac.ebi.pride.spectracluster.util.function.peak;

import uk.ac.ebi.pride.spectracluster.spectrum.IPeak;
import uk.ac.ebi.pride.spectracluster.spectrum.Peak;
import uk.ac.ebi.pride.spectracluster.util.function.IFunction;

import java.util.ArrayList;
import java.util.List;

/**
 * Transforms the peak intensities based on the formula
 * Inew = 1 + log(Iorg)
 * This function is taken from the spectral-archives algorithm.
 * Created by jg on 05.02.15.
 */
public class LogPeakIntensityTransformFunction implements IFunction<List<IPeak>, List<IPeak>> {
    @Override
    public List<IPeak> apply(List<IPeak> peaks) {
        List<IPeak> transformedIntensities = new ArrayList<IPeak>(peaks.size());

        for (IPeak originalPeak : peaks) {
            float transformedIntensity =  (originalPeak.getIntensity() != 0) ? 1 + (float) Math.log(originalPeak.getIntensity()) : 0;
            IPeak transformedPeak = new Peak(originalPeak.getMz(), transformedIntensity, originalPeak.getCount());
            transformedIntensities.add(transformedPeak);
        }

        return transformedIntensities;
    }
}

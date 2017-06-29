package uk.ac.ebi.pride.spectracluster.util.function.peak;

import uk.ac.ebi.pride.spectracluster.spectrum.IPeak;
import uk.ac.ebi.pride.spectracluster.spectrum.Peak;
import uk.ac.ebi.pride.spectracluster.util.function.IFunction;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

/**
 * This function replaces the peaks' intensities by their rank. The
 * highest peak gets the intensity of the number of peaks (if there are
 * 50 peaks, the highest peak gets intensity 50, the second highest 49
 * and so on). This is based on the method used in
 * Shao W, Zhu K, Lam H., Proteomics. 2013 Nov;13(22):3273-83
 *
 * Created by jg on 18.04.15.
 */
public class RankAsIntensityFunction implements IFunction<List<IPeak>, List<IPeak>> {
    @Override
    public List<IPeak> apply(List<IPeak> peaks) {
        List<Float> intensities = peaks.stream()
                .map(IPeak::getIntensity).sorted()
                .collect(Collectors.toCollection(() -> new ArrayList<>(peaks.size())));

        List<IPeak> changedPeaks = new ArrayList<>(peaks.size());

        for (IPeak p : peaks) {
            // get the rank
            int rank = intensities.indexOf(p.getIntensity()) + 1;

            changedPeaks.add(new Peak(p.getMz(), rank, p.getCount()));
        }

        return changedPeaks;
    }
}

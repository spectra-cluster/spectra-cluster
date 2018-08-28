package uk.ac.ebi.pride.spectracluster.normalizer;

import uk.ac.ebi.pride.spectracluster.spectrum.IPeak;
import uk.ac.ebi.pride.spectracluster.spectrum.Peak;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by jg on 24.01.15.
 */
public class SquaredSumIntensityNormalizer implements IIntensityNormalizer {
    public final String algorithmName = "SquaredSumIntensityNormalizer";
    public final String algorithmVersion = "1.0";

    @Override
    public double getNormalizedValue() {
        return 0; // TODO: this function is not very sensible
    }

    @Override
    public List<IPeak> normalizePeaks(List<IPeak> peaks) {
        // get the total intensity
        double totalSquaredIntensity = peaks.stream()
                .mapToDouble(p -> Math.abs(p.getIntensity()))
                .sum();

        double factor = Math.sqrt(totalSquaredIntensity);

        List<IPeak> normalizedPeaks = new ArrayList<>(peaks.size());

        for (IPeak p : peaks) {
            float normalizedIntensity = (float) ( Math.sqrt(p.getIntensity()) / factor );
            normalizedPeaks.add(new Peak(p.getMz(), normalizedIntensity, p.getCount()));
        }

        return normalizedPeaks;
    }

    @Override
    public String getName() {
        return algorithmName;
    }

    @Override
    public String getCurrentVersion() {
        return algorithmVersion;
    }
}

package uk.ac.ebi.pride.spectracluster.util.function.peak;

import uk.ac.ebi.pride.spectracluster.spectrum.IPeak;
import uk.ac.ebi.pride.spectracluster.util.comparator.PeakIntensityComparator;
import uk.ac.ebi.pride.spectracluster.util.function.IFunction;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;


/**
 * Filters the peaks based on a set total ion current. Keeps
 * the N highest peaks required to explain the set percentage
 * of the total spectrum intensity.
 * Created by jg on 04.02.15.
 */
public class FractionTICPeakFunction implements IFunction<List<IPeak>, List<IPeak>> {
    public final float DEFAULT_FRACTION_TOTAL_INTENSITY = 0.9F;
    /**
     * The total ion current that the N highest
     * peaks need to explain.
     */
    private float fractionTotalIntensity;

    public FractionTICPeakFunction() {
        this.fractionTotalIntensity = DEFAULT_FRACTION_TOTAL_INTENSITY;
    }

    public FractionTICPeakFunction(float percentageTotalIntensity) {
        this.fractionTotalIntensity = percentageTotalIntensity;
    }

    @Override
    public List<IPeak> apply(List<IPeak> peaks) {
        // calculate the total intensity
        double totalIntensity = 0;
        for (IPeak p : peaks)
            totalIntensity += p.getIntensity();

        // sort according to intensity (descending order)
        Collections.sort(peaks, PeakIntensityComparator.INSTANCE);

        List<IPeak> filteredPeaks = new ArrayList<IPeak>();
        double retainedIntensity = 0;

        for (IPeak p : peaks) {
            if (retainedIntensity / totalIntensity > fractionTotalIntensity)
                break;

            filteredPeaks.add(p);
            retainedIntensity += p.getIntensity();
        }

        return filteredPeaks;
    }

    public float getFractionTotalIntensity() {
        return fractionTotalIntensity;
    }
}

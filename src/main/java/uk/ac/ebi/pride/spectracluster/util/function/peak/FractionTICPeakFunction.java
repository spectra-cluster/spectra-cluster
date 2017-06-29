package uk.ac.ebi.pride.spectracluster.util.function.peak;

import uk.ac.ebi.pride.spectracluster.spectrum.IPeak;
import uk.ac.ebi.pride.spectracluster.util.comparator.PeakIntensityComparator;
import uk.ac.ebi.pride.spectracluster.util.function.IFunction;

import java.util.ArrayList;
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
    /**
     * Minimum number of peaks to keep
     */
    private int minimumNumberOfPeaks = 0;

    public FractionTICPeakFunction() {
        this.fractionTotalIntensity = DEFAULT_FRACTION_TOTAL_INTENSITY;
    }

    public FractionTICPeakFunction(float percentageTotalIntensity) {
        this.fractionTotalIntensity = percentageTotalIntensity;
    }

    public FractionTICPeakFunction(float fractionTotalIntensity, int minimumNumberOfPeaks) {
        this.fractionTotalIntensity = fractionTotalIntensity;
        this.minimumNumberOfPeaks = minimumNumberOfPeaks;
    }

    @Override
    public List<IPeak> apply(List<IPeak> peaks) {
        // calculate the total intensity
        double totalIntensity = peaks.stream().mapToDouble(IPeak::getIntensity).sum();

        // sort according to intensity (descending order)
        List<IPeak> sortedPeaks = new ArrayList<>(peaks);
        sortedPeaks.sort(PeakIntensityComparator.INSTANCE);

        List<IPeak> filteredPeaks = new ArrayList<>();
        double retainedIntensity = 0;
        int nPeaksRetained = 0;

        for (IPeak p : sortedPeaks) {
            if (retainedIntensity / totalIntensity > fractionTotalIntensity && nPeaksRetained >= minimumNumberOfPeaks)
                break;

            filteredPeaks.add(p);
            retainedIntensity += p.getIntensity();
            nPeaksRetained++;
        }

        return filteredPeaks;
    }

    public float getFractionTotalIntensity() {
        return fractionTotalIntensity;
    }

    public int getMinimumNumberOfPeaks() {
        return minimumNumberOfPeaks;
    }
}

package uk.ac.ebi.pride.spectracluster.similarity;

import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.stat.correlation.KendallsCorrelation;
import uk.ac.ebi.pride.spectracluster.spectrum.IPeak;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.util.comparator.PeakIntensityComparator;

import java.util.*;

/**
 * This SimilarityChecker assess the similarity between two spectra
 * by using the Kenall-Tau rank correlation coefficient of the intensities
 * of the matched peaks.
 * Created by jg on 23.02.15.
 */
public class IntensityRankCorrelation implements ISimilarityChecker {
    /**
     * The m/z tolerance to use in peak matching
     */
    private float peakMzTolerance = 0.5F;
    private boolean peakFiltering = true;

    private KendallsCorrelation kendallsCorrelation = new KendallsCorrelation();

    public IntensityRankCorrelation() {

    }

    public IntensityRankCorrelation(float peakMzTolerance) {
        this.peakMzTolerance = peakMzTolerance;
    }

    public IntensityRankCorrelation(float peakMzTolerance, boolean peakFiltering) {
        this.peakMzTolerance = peakMzTolerance;
        this.peakFiltering = peakFiltering;
    }

    @Override
    public double assessSimilarity(PeakMatches peakMatches) {
        // if there are no shared peaks, return 0 to indicate that it's random
        if (peakMatches.getNumberOfSharedPeaks() < 1)
            return 1;

        // only use the intensities
        double[] intensitiesSpec1 = extractPeakIntensities(peakMatches.getSharedPeaksSpec1());
        double[] intensitiesSpec2 = extractPeakIntensities(peakMatches.getSharedPeaksSpec2());

        double correlation = kendallsCorrelation.correlation(intensitiesSpec1, intensitiesSpec2);

        // convert correlation into probability using the distribution used in Pepitome
        // Normal Distribution with mean = 0 and SD^2 = 2(2k + 5)/9k(k − 1)
        double k = (double) peakMatches.getNumberOfSharedPeaks();
        // this cannot be calculated for only 1 shared peak
        if (k == 1)
            return 0;
        double sdSquare = (2 * (2 * k + 5)) / (9 * k * (k - 1) );
        double sd = Math.sqrt(sdSquare);
        NormalDistribution correlationDistribution = new NormalDistribution(0, sd);
        double probability = correlationDistribution.cumulativeProbability(correlation);

        return -Math.log(1 - probability);
    }

    @Override
    public double assessSimilarity(ISpectrum spectrum1, ISpectrum spectrum2) {
        // get similar peaks
        ISpectrum filteredSpectrum1, filteredSpectrum2;

        if (peakFiltering) {
            int nPeaks = calculateNPeaks(spectrum1.getPrecursorMz(), spectrum2.getPrecursorMz());
            if (nPeaks < 20)
                nPeaks = 20;

            filteredSpectrum1 = spectrum1.getHighestNPeaks(nPeaks);
            filteredSpectrum2 = spectrum2.getHighestNPeaks(nPeaks);
        }
        else {
            // simply disable filtering
            filteredSpectrum1 = spectrum1;
            filteredSpectrum2 = spectrum2;
        }

        PeakMatches peakMatches = SimilarityUtilities.getSharedPeaksAsMatches(filteredSpectrum1, filteredSpectrum2, peakMzTolerance);

        return assessSimilarity(peakMatches);
    }

    private double[] extractPeakIntensities(IPeak[] peaks) {
        double[] intensities = new double[peaks.length];

        for (int i = 0; i < peaks.length; i++) {
            intensities[i] = (double) peaks[i].getIntensity();
        }

        return intensities;
    }

    @Override
    public boolean isPeakFiltering() {
        return peakFiltering;
    }

    @Override
    public void setPeakFiltering(boolean peakFiltering) {
        this.peakFiltering = peakFiltering;
    }

    @Override
    public String getName() {
        return null;
    }

    @Override
    public String getCurrentVersion() {
        return null;
    }

    /**
     * Calculate k by using the precursor m/z / 50.
     *
     * @param precursor1
     * @param precursor2
     * @return
     */
    private int calculateNPeaks(Float precursor1, Float precursor2) {
        // if any of the required values is missing, return 15
        if (precursor1 == null || precursor2 == null)
            return 15;

        // use m/z / 50

        return (int) ((precursor1 / 50 + precursor2 / 50) / 2);
    }
}

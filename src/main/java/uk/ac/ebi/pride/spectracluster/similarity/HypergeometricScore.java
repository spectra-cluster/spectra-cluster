package uk.ac.ebi.pride.spectracluster.similarity;

import org.apache.commons.math3.distribution.HypergeometricDistribution;
import uk.ac.ebi.pride.spectracluster.spectrum.IPeak;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.util.Defaults;

import java.util.List;

/**
 * This is an implementation of the Spectral Comparison
 * function used in Pepitope. It is based on the hypergeometric
 * probability that the observed similar m/z values are a
 * random match.
 * <p/>
 * Created by jg on 15.01.15.
 */
public class HypergeometricScore implements ISimilarityChecker {
    public static final String algorithmName = "Hypergeometric Exact Test";
    public static final String algorithmVersion = "0.1";

    public static final boolean DEFAULT_PEAK_FILTERING = true;

    private boolean peakFiltering;

    /**
     * The tolerance in m/z units used to match peaks
     */
    protected float fragmentIonTolerance;

    public HypergeometricScore() {
        this(Defaults.getFragmentIonTolerance(), DEFAULT_PEAK_FILTERING);
    }

    public HypergeometricScore(float fragmentIonTolerance) {
        this.fragmentIonTolerance = fragmentIonTolerance;
        this.peakFiltering = DEFAULT_PEAK_FILTERING;
    }

    public HypergeometricScore(float fragmentIonTolerance, boolean peakFiltering) {
        this.fragmentIonTolerance = fragmentIonTolerance;
        this.peakFiltering = peakFiltering;
    }

    @Override
    public double assessSimilarity(IPeakMatches peakMatches) {
        // if there are no shared peaks, return 0 to indicate that it's random
        if (peakMatches.getNumberOfSharedPeaks() < 1)
            return 1;

        int numberOfBins = calculateNumberOfBins(peakMatches);

        return calculateSimilarityScore(peakMatches.getNumberOfSharedPeaks(),
                peakMatches.getSpectrumOne().getPeaksCount(),
                peakMatches.getSpectrumTwo().getPeaksCount(),
                numberOfBins);
    }

    public double assessSimilarityAsPValue(IPeakMatches peakMatches) {
        // if there are no shared peaks, return 0 to indicate that it's random
        if (peakMatches.getNumberOfSharedPeaks() < 1)
            return 1;

        int numberOfBins = calculateNumberOfBins(peakMatches);

        return calculateSimilarityProbablity(peakMatches.getNumberOfSharedPeaks(),
                peakMatches.getSpectrumOne().getPeaksCount(),
                peakMatches.getSpectrumTwo().getPeaksCount(),
                numberOfBins);
    }

    protected int calculateNumberOfBins(IPeakMatches peakMatches) {
        List<IPeak> peaks1 = peakMatches.getSpectrumOne().getPeaks();
        List<IPeak> peaks2 = peakMatches.getSpectrumTwo().getPeaks();

        // set the maximum shared m/z value
        float minMz, maxMz; // minimum and maximum overlapping m/z

        if (peaks1.get(0).getMz() < peaks2.get(0).getMz()) {
            minMz = peaks1.get(0).getMz();
        } else {
            minMz = peaks2.get(0).getMz();
        }

        if (peaks1.get(peaks1.size() - 1).getMz() > peaks2.get(peaks2.size() - 1).getMz()) {
            maxMz = peaks1.get(peaks1.size() - 1).getMz();
        } else {
            maxMz = peaks2.get(peaks2.size() - 1).getMz();
        }

        int numberOfBins = Math.round((maxMz - minMz) / fragmentIonTolerance);

        // cannot be assessed
        if (numberOfBins < 1) {
            return 0;
        }

        if (numberOfBins < peaks1.size()|| numberOfBins < peaks2.size()) {
            return 0;
        }

        return numberOfBins;
    }

    protected double calculateSimilarityProbablity(int numberOfSharedPeaks, int numberOfPeaksFromSpec1, int numberOfPeaksFromSpec2, int numberOfBins) {
        if (numberOfBins < 1) {
            return 1;
        }

        // ToDo: @jgriss In peptidome manuscript, the number of successes and the sample size are the same, was it a mistake from them?
        HypergeometricDistribution hypergeometricDistribution = new HypergeometricDistribution(numberOfBins, numberOfPeaksFromSpec1, numberOfPeaksFromSpec2);

        double hgtScore = 0; // summed probability of finding more peaks
        for (int nFoundPeaks = numberOfSharedPeaks + 1; nFoundPeaks <= numberOfPeaksFromSpec2; nFoundPeaks++) {
            hgtScore += hypergeometricDistribution.probability(nFoundPeaks);
        }

        if (hgtScore == 0) {
            return 1;
        }

        return hgtScore;
    }

    protected double calculateSimilarityScore(int numberOfSharedPeaks, int numberOfPeaksFromSpec1, int numberOfPeaksFromSpec2, int numberOfBins) {
        double pValue = calculateSimilarityProbablity(numberOfSharedPeaks, numberOfPeaksFromSpec1, numberOfPeaksFromSpec2, numberOfBins);

        return -Math.log(pValue);
    }

    @Override
    public double assessSimilarity(ISpectrum spectrum1, ISpectrum spectrum2) {
        IPeakMatches peakMatches = PeakMatchesUtilities.getSharedPeaksAsMatches(spectrum1, spectrum2, fragmentIonTolerance, peakFiltering);
        return assessSimilarity(peakMatches);
    }

    @Override
    public String getName() {
        return algorithmName;
    }

    @Override
    public String getCurrentVersion() {
        return algorithmVersion;
    }

    @Override
    public void setFragmentIonTolerance(float fragmentIonTolerance) {
        this.fragmentIonTolerance = fragmentIonTolerance;
    }

    @Override
    public float getFragmentIonTolerance() {
        return fragmentIonTolerance;
    }

    @Override
    public boolean isPeakFiltering() {
        return peakFiltering;
    }

    @Override
    public void setPeakFiltering(boolean peakFiltering) {
        this.peakFiltering = peakFiltering;
    }
}

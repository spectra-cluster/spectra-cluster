package uk.ac.ebi.pride.spectracluster.similarity;

import org.apache.commons.math3.distribution.HypergeometricDistribution;
import uk.ac.ebi.pride.spectracluster.spectrum.IPeak;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;

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

    public static final float DEFAULT_PEAK_MZ_TOLERANCE = 0.5F;
    public static final boolean DEFAULT_PEAK_FILTERING = true;

    private boolean peakFiltering;

    /**
     * The tolerance in m/z units used to match peaks
     */
    protected float peakMzTolerance;

    public HypergeometricScore() {
        this(DEFAULT_PEAK_MZ_TOLERANCE, DEFAULT_PEAK_FILTERING);
    }

    public HypergeometricScore(float peakMzTolerance) {
        this.peakMzTolerance = peakMzTolerance;
    }

    public HypergeometricScore(float peakMzTolerance, boolean peakFiltering) {
        this.peakMzTolerance = peakMzTolerance;
        this.peakFiltering = peakFiltering;
    }

    @Override
    public double assessSimilarity(IPeakMatches peakMatches) {
        // if there are no shared peaks, return 0 to indicate that it's random
        if (peakMatches.getNumberOfSharedPeaks() < 1)
            return 1;

        // set the maximum shared m/z value
        float minMz, maxMz; // minimum and maximum overlapping m/z

        List<IPeak> peaks1 = peakMatches.getSpectrumOne().getPeaks();
        IPeak[] peaksSpec1 = peaks1.toArray(new IPeak[peaks1.size()]);

        List<IPeak> peaks2 = peakMatches.getSpectrumTwo().getPeaks();
        IPeak[] peaksSpec2 = peaks2.toArray(new IPeak[peaks2.size()]);

        if (peaksSpec1[0].getMz() < peaksSpec2[0].getMz()) {
            minMz = peaksSpec2[0].getMz();
        } else {
            minMz = peaksSpec1[0].getMz();
        }

        if (peaksSpec1[peaksSpec1.length - 1].getMz() > peaksSpec2[peaksSpec2.length - 1].getMz()) {
            maxMz = peaksSpec2[peaksSpec2.length - 1].getMz();
        } else {
            maxMz = peaksSpec1[peaksSpec1.length - 1].getMz();
        }

        int numberOfBins = Math.round((maxMz - minMz) / peakMzTolerance);

        // cannot be assessed
        if (numberOfBins < 1) {
            return 0;
        }

        if (numberOfBins < peaksSpec1.length || numberOfBins < peaksSpec2.length) {
            return 0;
        }

        return calculateSimilarityScore(peakMatches.getNumberOfSharedPeaks(), peaksSpec1.length, peaksSpec2.length, numberOfBins);


    }

    protected double calculateSimilarityScore(int numberOfSharedPeaks, int numberOfPeaksFromSpec1, int numberOfPeaksFromSpec2, int numberOfBins) {
        // ToDo: @jgriss In peptidome manuscript, the number of successes and the sample size are the same, was it a mistake from them?
        HypergeometricDistribution hypergeometricDistribution = new HypergeometricDistribution(numberOfBins, numberOfPeaksFromSpec1, numberOfPeaksFromSpec2);

        double hgtScore = 0; // summed probability of finding more peaks
        for (int nFoundPeaks = numberOfSharedPeaks + 1; nFoundPeaks <= numberOfPeaksFromSpec2; nFoundPeaks++) {
            hgtScore += hypergeometricDistribution.probability(nFoundPeaks);
        }

        if (hgtScore == 0) {
            return 0;
        }

        return -Math.log(hgtScore);
    }

    @Override
    public double assessSimilarity(ISpectrum spectrum1, ISpectrum spectrum2) {
        ISpectrum filteredSpectrum1, filteredSpectrum2;

        if (peakFiltering) {
            int nPeaks = PeakMatchesUtilities.calculateNPeaks(spectrum1.getPrecursorMz(), spectrum2.getPrecursorMz());
            if (nPeaks < 20)
                nPeaks = 20;

            filteredSpectrum1 = spectrum1.getHighestNPeaks(nPeaks);
            filteredSpectrum2 = spectrum2.getHighestNPeaks(nPeaks);
        } else {
            filteredSpectrum1 = spectrum1;
            filteredSpectrum2 = spectrum2;
        }

        IPeakMatches peakMatches = PeakMatchesUtilities.getSharedPeaksAsMatches(filteredSpectrum1, filteredSpectrum2, peakMzTolerance);

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

    public float getPeakMzTolerance() {
        return peakMzTolerance;
    }

    public void setPeakMzTolerance(float peakMzTolerance) {
        this.peakMzTolerance = peakMzTolerance;
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

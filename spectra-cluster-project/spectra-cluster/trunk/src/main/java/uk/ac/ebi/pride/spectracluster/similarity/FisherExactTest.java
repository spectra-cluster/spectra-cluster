package uk.ac.ebi.pride.spectracluster.similarity;

import org.apache.commons.math3.distribution.HypergeometricDistribution;
import uk.ac.ebi.pride.spectracluster.spectrum.IPeak;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;

import java.util.List;

/**
 * This SimilarityChecker is based on the hypergeometric
 * probability that the observed similar m/z values are a
 * random match. It only calculates a point probability.
 *
 * Created by jg on 15.01.15.
 */
public class FisherExactTest implements ISimilarityChecker {
    public static final String algorithmName = "Fisher Exact Test";
    public static final String algorithmVersion = "0.1";
    private boolean peakFiltering = true;

    /**
     * The tolerance in m/z units used to match peaks
     */
    protected float peakMzTolerance = 0.5F;

    public FisherExactTest() {

    }

    public FisherExactTest(float peakMzTolerance) {
        this.peakMzTolerance = peakMzTolerance;
    }

    public FisherExactTest(float peakMzTolerance, boolean peakFiltering) {
        this.peakFiltering = peakFiltering;
        this.peakMzTolerance = peakMzTolerance;
    }

    @Override
    public double assessSimilarity(PeakMatches peakMatches) {
        // if there are no shared peaks, return 0 to indicate that it's random
        if (peakMatches.getNumberOfSharedPeaks() < 1)
            return 1;

        // set the maximum shared m/z value
        float minMz, maxMz; // minimum and maximum overlapping m/z

        IPeak[] peaksSpec1 = peakMatches.getSpectrum1().getPeaks().toArray(new IPeak[peakMatches.getSpectrum1().getPeaks().size()]);
        IPeak[] peaksSpec2 = peakMatches.getSpectrum2().getPeaks().toArray(new IPeak[peakMatches.getSpectrum2().getPeaks().size()]);

        if (peaksSpec1[0].getMz() < peaksSpec2[0].getMz()) {
            minMz = peaksSpec2[0].getMz();
        }
        else {
            minMz = peaksSpec1[0].getMz();
        }

        if (peaksSpec1[peaksSpec1.length - 1].getMz() > peaksSpec2[peaksSpec2.length - 1].getMz()) {
            maxMz = peaksSpec2[peaksSpec2.length - 1].getMz();
        }
        else {
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

        HypergeometricDistribution hypergeometricDistribution = new HypergeometricDistribution(
                numberOfBins, peaksSpec1.length, peaksSpec2.length);

        // point probability
        double hgtScore = hypergeometricDistribution.probability(peakMatches.getNumberOfSharedPeaks());

        return -Math.log(hgtScore);
    }

    @Override
    public double assessSimilarity(ISpectrum spectrum1, ISpectrum spectrum2) {
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

    public boolean isPeakFiltering() {
        return peakFiltering;
    }

    public void setPeakFiltering(boolean peakFiltering) {
        this.peakFiltering = peakFiltering;
    }
}

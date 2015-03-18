package uk.ac.ebi.pride.spectracluster.similarity;

import org.apache.commons.math3.distribution.HypergeometricDistribution;
import uk.ac.ebi.pride.spectracluster.spectrum.IPeak;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;

import java.util.List;

/**
 * This version of the HypergeometricScore uses the total number of
 * bins from spectrum 1 as a population size. The original version uses
 * the total number of bins within the overlapping region of both spectra.
 * @author Johannes Griss
 * @author Rui Wang
 */
public class HypergeometricScoreDiffPopSize extends HypergeometricScore {
    public HypergeometricScoreDiffPopSize() {
    }

    public HypergeometricScoreDiffPopSize(float peakMzTolerance) {
        super(peakMzTolerance);
    }

    public HypergeometricScoreDiffPopSize(float peakMzTolerance, boolean peakFiltering) {
        super(peakMzTolerance, peakFiltering);
    }

    @Override
    public double assessSimilarity(IPeakMatches peakMatches) {
        // if there are no shared peaks, return 0 to indicate that it's random
        if (peakMatches.getNumberOfSharedPeaks() < 1)
            return 1;

        List<IPeak> peaksSpectrumOne = peakMatches.getSpectrumOne().getPeaks();

        float minMz = peaksSpectrumOne.get(0).getMz();
        float maxMz = peaksSpectrumOne.get(peaksSpectrumOne.size() - 1).getMz();

        int numberOfBins = Math.round((maxMz - minMz) / peakMzTolerance);

        return calculateSimilarityScore(peakMatches.getNumberOfSharedPeaks(),
                peakMatches.getSpectrumOne().getPeaksCount(),
                peakMatches.getSpectrumTwo().getPeaksCount(),
                numberOfBins);
    }
}

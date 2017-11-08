package uk.ac.ebi.pride.spectracluster.similarity;

import cern.jet.random.Normal;
import cern.jet.random.engine.RandomEngine;
import org.apache.commons.math3.stat.correlation.KendallsCorrelation;
import uk.ac.ebi.pride.spectracluster.spectrum.IPeak;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.util.Defaults;

import java.util.List;

/**
 * This SimilarityChecker assess the similarity between two spectra
 * by using the Kenall-Tau rank correlation coefficient of the intensities
 * of the matched peaks.
 * Created by jg on 23.02.15.
 *
 * @author Yasset Perez-Riverol
 */
public class IntensityRankCorrelation implements ISimilarityChecker {
    public final static boolean DEFAULT_PEAK_FILTERING = false;
    protected final RandomEngine randomEngine = RandomEngine.makeDefault();

    /**
     * The m/z tolerance to use in peak matching
     */
    protected float fragmentIonTolerance;
    protected boolean peakFiltering;

    private SerializableKendallsCorrelation kendallsCorrelation = new SerializableKendallsCorrelation();

    public IntensityRankCorrelation() {
        this(Defaults.getFragmentIonTolerance(), DEFAULT_PEAK_FILTERING);
    }

    public IntensityRankCorrelation(float fragmentIonTolerance) {
        this(fragmentIonTolerance, DEFAULT_PEAK_FILTERING);
    }

    public IntensityRankCorrelation(float fragmentIonTolerance, boolean peakFiltering) {
        this.fragmentIonTolerance = fragmentIonTolerance;
        this.peakFiltering = peakFiltering;
    }

    /**
     * This function provides a p-value for the similarity between two different spectrum
     * based on kendallsCorrelation
     *
     * @param peakMatches Number of peaks that matches between two spectra
     * @return p-value of the similarity
     */
    public double assessSimilarityAsPValue(IPeakMatches peakMatches) {
        // if there are no shared peaks, return 1 to indicate that it's random
        if (peakMatches.getNumberOfSharedPeaks() < 1)
            return 1;

        // only use the intensities
        double[] intensitiesSpec1 = extractPeakIntensities(peakMatches.getSharedPeaksFromSpectrumOne());
        double[] intensitiesSpec2 = extractPeakIntensities(peakMatches.getSharedPeaksFromSpectrumTwo());

        double correlation = kendallsCorrelation.correlation(intensitiesSpec1, intensitiesSpec2);

        // if the correlation cannot be calculated, assume that there is none
        if (Double.isNaN(correlation)) {
            return 1;
        }

        // convert correlation into probability using the distribution used in Peptidome
        // Normal Distribution with mean = 0 and SD^2 = 2(2k + 5)/9k(k − 1)
        double k = (double) peakMatches.getNumberOfSharedPeaks();

        // this cannot be calculated for only 1 shared peak
        if (k == 1)
            return 1;

        double sdSquare = (2 * (2 * k + 5)) / (9 * k * (k - 1) );
        double sd = Math.sqrt(sdSquare);

        Normal normal = new Normal(0, sd, randomEngine);
        double probability = normal.cdf(correlation);

        return 1 - probability;
    }

    @Override
    public double assessSimilarity(IPeakMatches peakMatches) {
        double pValue = assessSimilarityAsPValue(peakMatches);
        return -Math.log(pValue);
    }

    @Override
    public double assessSimilarity(ISpectrum spectrum1, ISpectrum spectrum2) {
        IPeakMatches peakMatches = PeakMatchesUtilities.getSharedPeaksAsMatches(spectrum1, spectrum2, fragmentIonTolerance, peakFiltering);
        return assessSimilarity(peakMatches);
    }

    private double[] extractPeakIntensities(List<IPeak> peaks) {
        double[] intensities = new double[peaks.size()];

        for (int i = 0; i < peaks.size(); i++) {
            intensities[i] = (double) peaks.get(i).getIntensity();
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

    @Override
    public void setFragmentIonTolerance(float fragmentIonTolerance) {
        this.fragmentIonTolerance = fragmentIonTolerance;
    }

    @Override
    public float getFragmentIonTolerance() {
        return fragmentIonTolerance;
    }
}

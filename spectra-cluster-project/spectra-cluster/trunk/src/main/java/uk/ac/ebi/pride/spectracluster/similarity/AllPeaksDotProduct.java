package uk.ac.ebi.pride.spectracluster.similarity;


import uk.ac.ebi.pride.spectracluster.spectrum.IPeak;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.spectrum.Peak;

import java.util.ArrayList;
import java.util.List;


/**
 * Assesses the similarity between two
 * spectra using the normalized
 * dot-product as implemented by Frank et al. (2008) JPR
 * <p/>
 * This implementation uses all peaks of the two spectra instead of
 * only the k-highest ones.
 * <p/>
 * uk.ac.ebi.pride.spectracluster.similarity.FrankEtAlDotProduct
 *
 * @author jg
 *         <p/>
 */
// this seems to be for testing only SLewis
@Deprecated
public class AllPeaksDotProduct implements ISimilarityChecker {
    private final static String version = "0.1";


    // add a nonsense park
    private static final IPeak LAST_PEAK = new Peak(Float.MAX_VALUE, 0);

    private double similarityMZRange;

    public AllPeaksDotProduct(double similarityMZRange) {
        this.similarityMZRange = similarityMZRange;
    }

    /**
     * return a name which should not change
     *
     * @return !null name
     */
    @Override
    public String getName() {
        return getClass().getSimpleName();
    }

    /**
     * return a version number - this may be updated over time
     *
     * @return !null version
     */
    @Override
    public String getCurrentVersion() {
        return version;
    }

    /**
     * Assesses the spectra's similarity using
     * the normalized dot-product
     */
    @Override
    public double assessSimilarity(ISpectrum spectrum1, ISpectrum spectrum2) {

        // get the k highest peaks1 from every spectrum
        double sumSquareIntensity1 = spectrum1.getSumSquareIntensity();

        // the collection is immutable we need to build a new one
        List<IPeak> peaksSpec1 = new ArrayList<IPeak>(spectrum1.getPeaks());
        peaksSpec1.add(LAST_PEAK); //add a peak we will not use
        IPeak[] peaks1 = peaksSpec1.toArray(new IPeak[peaksSpec1.size()]);

        double sumSquareIntensity2 = spectrum2.getSumSquareIntensity();
        // the collection is immutable we need to build a new one
        List<IPeak> peaksSpec2 = new ArrayList<IPeak>(spectrum2.getPeaks());
        peaksSpec2.add(LAST_PEAK); //add a peak we will not use
        IPeak[] peaks2 = peaksSpec2.toArray(new IPeak[peaksSpec2.size()]);

        double mzRange = similarityMZRange;
        boolean lastIsT = false;
        int t = 0;
        int e = 0;
        double dotProduct = 0.0;

        while (t < peaks1.length - 1 && e < peaks2.length - 1) {
            IPeak peak1 = peaks1[t];
            double mz1 = peak1.getMz();
            IPeak peak2 = peaks2[e];
            double mz2 = peak2.getMz();

            double mass_difference = mz2 - mz1;

            if (Math.abs(mass_difference) <= mzRange) {
                // also calculate the difference for the next t and e peaks
                IPeak nextPeak1 = peaks1[t + 1];
                IPeak nextPeak2 = peaks2[e + 1];
                double mass_difference_nextT = mz2 - nextPeak1.getMz();
                double mass_difference_nextE = nextPeak2.getMz() - mz1;
                double mass_difference_nextTE = nextPeak2.getMz() - nextPeak1.getMz();

                // use the next spectrum in E if it's a better match
                if (Math.abs(mass_difference_nextE) < Math.abs(mass_difference) &&
                        Math.abs(mass_difference_nextE) < Math.abs(mass_difference_nextT) &&
                        Math.abs(mass_difference_nextE) < Math.abs(mass_difference_nextTE)) {
                    e++;
                    peak2 = nextPeak2;
                    mass_difference = mass_difference_nextE;
                }

                // use the next spectrum in T if it's a better match
                if (Math.abs(mass_difference_nextT) < Math.abs(mass_difference) &&
                        Math.abs(mass_difference_nextT) < Math.abs(mass_difference_nextE) &&
                        Math.abs(mass_difference_nextT) < Math.abs(mass_difference_nextTE)) {
                    t++;
                    peak1 = nextPeak1;
                    mass_difference = mass_difference_nextT;
                }

                // do the matching
                double match1 = convertIntensity(peak1);
                double match2 = convertIntensity(peak2);

                dotProduct += match1 * match2;

                // increment both counters since both peaks must not be compared again
                e++;
                t++;
                continue;
            }

            if (mass_difference == 0) {
                if (lastIsT) {
                    e++;
                    lastIsT = false;
                } else {
                    t++;
                    lastIsT = true;
                }
            } else {
                if (mass_difference < 0) {
                    e++;
                } else {
                    t++;
                }
            }

        }
        // normalize the dot product
        double denom = Math.sqrt(sumSquareIntensity1 * sumSquareIntensity2);
        if (denom == 0)
            return 0;
        double normalizedDotProduct = dotProduct / denom;

        if (normalizedDotProduct > 1)
            throw new IllegalStateException("Dot product > 1. This is mathematically not possible.");

        return normalizedDotProduct;
    }

    /**
     * who knows why Johannes does this but we can as well
     * todo @rw: double check this wit Johannes
     */
    private double convertIntensity(IPeak p1) {
        double intensity = p1.getIntensity();
        if (intensity == 0)
            return 0;
        return 1 + Math.log(intensity);
    }
}

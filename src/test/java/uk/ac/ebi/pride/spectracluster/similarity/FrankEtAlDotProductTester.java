package uk.ac.ebi.pride.spectracluster.similarity;


import uk.ac.ebi.pride.spectracluster.spectrum.IPeak;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;

/**
 * Assesses the similarity between two
 * spectra using the normalized
 * dot-product as implemented by Frank et al. (2008) JPR
 * This implementation uses only the k highest peaks.
 * k is calculated by taking 15 peaks per 1000 Da
 * peptide mass. Furthermore, the vectors for the
 * dot-product are filled with the 1+ln(I) where I
 * is the peak's normalized intensity.
 *
 * @author jg
 *         todo: this class needs to be reviewed
 */
public class FrankEtAlDotProductTester implements ISimilarityChecker {

    public static final double DEFAULT_SIMILARITY_THRESHOLD = 0.6;

    /**
     * The versions available from this algorithm. The only
     * difference is the way k is calculated.
     *
     * @author jg
     */
    public enum AlgorithmVersion {
        JPR_2008, NAT_METH_2011
    }

    private double mzRange = 0.5;

    private boolean useClosestPeak = true;

    /**
     * The algorithm version to use. By
     * default the version described in
     * Nature Methods 2011 will be used.
     */
    private AlgorithmVersion version = AlgorithmVersion.NAT_METH_2011;

    /**
     * Use Defaults which builds with reflection
     * Set the class with Defaults.setSimilarityCheckerClass
     */
    public FrankEtAlDotProductTester() {
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
        return version.toString();
    }

    @Override
    public double assessSimilarity(IPeakMatches peakMatches) {
        throw new UnsupportedOperationException();
    }

    /**
     * Assesses the spectra's similarity using
     * the normalized dot-product
     */
    @Override
    public double assessSimilarity(ISpectrum spectrum1, ISpectrum spectrum2) {
        // initialize the number of peaks to use with 15
        int k = 15;
        switch (version) {
            case NAT_METH_2011:
                k = calculateK2011(spectrum1.getPrecursorMz(), spectrum2.getPrecursorMz());
                break;
            case JPR_2008:
                k = calculateK2008(spectrum1.getPrecursorMz(), spectrum2.getPrecursorMz(), spectrum1.getPrecursorCharge(), spectrum2.getPrecursorCharge());
                break;
        }

        // get the k highest peaks from every spectrum
        //List<IPeak> kHighestPeaks1 = getHighestPeaks(new ArrayList<IPeak>(spectrum1.getPeaks()), k);
        //List<IPeak> kHighestPeaks2 = getHighestPeaks(new ArrayList<IPeak>(spectrum2.getPeaks()), k);
        List<IPeak> kHighestPeaks1 = spectrum1.getHighestNPeaks(k).getPeaks();
        List<IPeak> kHighestPeaks2 = spectrum2.getHighestNPeaks(k).getPeaks();

        // order the two peak lists based on their m/z values
        //Collections.sort(kHighestPeaks1, PeakMzComparator.getInstance());
        //Collections.sort(kHighestPeaks2, PeakMzComparator.getInstance());

        // create two intensity vectors
        List<Double> intensities1 = new ArrayList<Double>(k * 2);
        List<Double> intensities2 = new ArrayList<Double>(k * 2);

        // indicates the last item in the k2HighestPeakList that was merged
        int lastIndex2 = 0;
        HashSet<Integer> addedPeaks2 = new HashSet<Integer>();

        for (IPeak peak1 : kHighestPeaks1) {
            // add the intensity to the intensity array of spectrum 1
            double intensity = peak1.getIntensity();
            double intensity2 = 1 + Math.log(intensity);
            intensities1.add(intensity2);

            double mz1 = peak1.getMz();

            // get the indexes of the comparable masses from peak list 2
            List<Integer> comparableIndexes = new ArrayList<Integer>(3);

            // TODO: @Steve One thing that's currently missing in this test is that it only checks
            // the peaks from spectrum 2 if there's a better fit. Ideally, the code should
            // only check spectrum 1 if there's a better fitting peak. For today this
            // code's enough as a proof of principle. In the final version I'd want to
            // check both spectra.

            // TODO: @Steve Here peaks within the tolerance are selected.
            for (int i = lastIndex2; i < kHighestPeaks2.size(); i++) {
                // make sure the peak wasn't used before
                if (addedPeaks2.contains(i))
                    continue;

                // compare the m/z
                double mz2 = kHighestPeaks2.get(i).getMz();

                // check if they are unique
                if (mz2 >= mz1 - mzRange && mz2 <= mz1 + mzRange) {
                    comparableIndexes.add(i);
                }

                // make sure the m/z isn't too big already
                if (mz2 > mz1 + mzRange)
                    break;
            }

            // get the comparable mass closest to the current one
            int closestIndex = -1;
            double closestDiff = 100;

            // TODO: @Steve here the closest matching peak is selected
            for (Integer i : comparableIndexes) {
                IPeak iPeak = kHighestPeaks2.get(i);
                double mz2 = iPeak.getMz();
                double diff = Math.abs(mz1 - mz2);

                if (diff < closestDiff) {
                    closestIndex = i;
                    closestDiff = diff;
                    // only use the first peak if not checking for closest peaks
                    if (!useClosestPeak)
                        break;
                }
            }

            // set the intensity 2
            if (closestIndex >= 0) {
                IPeak iPeak = kHighestPeaks2.get(closestIndex);
                double intensity1 = iPeak.getIntensity();
                intensities2.add(1 + Math.log(intensity1));
                addedPeaks2.add(closestIndex);
            } else {
                intensities2.add(0.0);
            }
        }

        // add the intensities for the second spectrum
        for (int i = 0; i < kHighestPeaks2.size(); i++) {
            if (addedPeaks2.contains(i))
                continue;

            // the peak doesn't exist in spectrum 1 as this was already checked
            intensities1.add(0.0);

            IPeak peak2 = kHighestPeaks2.get(i);

//                if(debugOutput != null)   {
//                       String fmt = String.format("%8.3f %8.3f  %8.3f %8.3f \n", peak1.getMz(), match1, peak2.getMz(), match2);
//                     try {
//                         debugOutput.append(fmt);
//                     }
//                     catch (IOException e1) {
//                         throw new RuntimeException(e1);
//
//                     }
//                  }

            intensities2.add(1 + Math.log(peak2.getIntensity()));
        }

        // make sure both intensities have the same size
        if (intensities1.size() != intensities2.size())
            throw new IllegalStateException("Different sizes of intensity arrays encountered.");

        // calculate the dot product
        double dotProduct = 0;
        double sumSquareIntensity1 = 0;
        double sumSquareIntensity2 = 0;
        int numberMatches = 0;
        for (int i = 0; i < intensities1.size(); i++) {
            double i1 = intensities1.get(i);
            double i2 = intensities2.get(i);
            //         String fmt  = String.format("%8.3f %8.3f", i1, i2);
            //        System.out.println(fmt);
            dotProduct += i1 * i2;
            numberMatches++;
            sumSquareIntensity1 += Math.pow(i1, 2);
            sumSquareIntensity2 += Math.pow(i2, 2);
        }

        // normalize the dot product
        double denom = Math.sqrt(sumSquareIntensity1 * sumSquareIntensity2);
        if (denom == 0)
            return 0;
        double normalizedDotProduct = dotProduct / denom;

        //     System.out.println("Old Spectrum matched " + numberMatches );
        return normalizedDotProduct;
    }

    /**
     * Calculate k by using 15 per 1000 Da of
     * peptide mass.
     *
     * @param precursor1
     * @param precursor2
     * @param charge1
     * @param charge2
     * @return
     */
    private int calculateK2008(Float precursor1, Float precursor2,
                               Integer charge1, Integer charge2) {
        // if any of the required values is missing, return 15
        if (precursor1 == null || precursor2 == null || charge1 == null || charge2 == null || charge1 <= 0 || charge2 <= 0)
            return 15;

        // take 15 peaks / 1000Da peptide mass
        double peptideMass = (precursor1 * charge1 + precursor2 * charge2) / 2;

        int k = 15 * (int) (peptideMass / 1000);

        if (peptideMass % 1000 > 0)
            k += 15;

        return k;
    }

    /**
     * Calculate k by using the precursor m/z / 50.
     *
     * @param precursor1
     * @param precursor2
     * @return
     */
    private int calculateK2011(Float precursor1, Float precursor2) {
        // if any of the required values is missing, return 15
        if (precursor1 == null || precursor2 == null)
            return 15;

        // use m/z / 50
        int k = (int) ((precursor1 / 50 + precursor2 / 50) / 2);

        return k;
    }

    public void setMzRange(float mzRange) {
        this.mzRange = mzRange;
    }

    public void setVersion(AlgorithmVersion version) {
        this.version = version;
    }

    public boolean isUseClosestPeak() {
        return useClosestPeak;
    }

    public void setUseClosestPeak(boolean useClosestPeak) {
        this.useClosestPeak = useClosestPeak;
    }

    @Override
    public boolean isPeakFiltering() {
        return true;
    }

    @Override
    public void setPeakFiltering(boolean peakFiltering) {

    }

    @Override
    public void setFragmentIonTolerance(float fragmentIonTolerance) {

    }

    @Override
    public float getFragmentIonTolerance() {
        return 0;
    }
}
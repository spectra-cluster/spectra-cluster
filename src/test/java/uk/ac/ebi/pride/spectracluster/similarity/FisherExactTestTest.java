package uk.ac.ebi.pride.spectracluster.similarity;

import org.junit.Before;
import org.junit.Test;
import uk.ac.ebi.pride.spectracluster.io.ParserUtilities;
import uk.ac.ebi.pride.spectracluster.normalizer.SquaredSumIntensityNormalizer;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.spectrum.Spectrum;

import java.io.InputStreamReader;
import java.io.LineNumberReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Created by jg on 17.01.15.
 */
public class FisherExactTestTest {
    public final String testFilename = "spectra_400.0_4.0.mgf";
    //public final String testFilename = "PBMC_3_stim_dex_sn.mgf";

    //public final String[] testFilenames = {"PBMC_1_stim_dex_sn.mgf.identified", "PBMC_2_stim_dex_sn.mgf.identified", "PBMC_3_stim_dex_sn.mgf.identified"};
    public final String[] testFilenames = {testFilename};

    List<ISpectrum> testSpectra;
    ISpectrum[] normalizedTestSpectra;

    @Before
    public void setup() throws Exception {
        testSpectra = new ArrayList<>();

        for (String filename : testFilenames) {
            ISpectrum[] spectra = ParserUtilities.readMGFScans(new LineNumberReader(new InputStreamReader(FisherExactTestTest.class.getClassLoader().getResourceAsStream(filename))));
            testSpectra.addAll(Arrays.asList(spectra));
        }

        // set the idenfied sequences
        if (testFilename.startsWith("PBMC_")) {
            for (ISpectrum s : testSpectra) {
                String id = s.getId();
                int index = id.indexOf("splib_sequence=");
                int endIndex = id.indexOf(",", index);

                if (index < 0 || endIndex < 0) {
                    throw new Exception("Failed to extract identification for " + id);
                }

                String sequence = id.substring(index + 15, endIndex);
                s.setProperty("splib_sequence", sequence);

                index = id.indexOf("peptideR2=");
                endIndex = id.indexOf(",", index);
                sequence = id.substring(index + 10, endIndex);
                s.setProperty("sequenceR2", sequence);
            }
        }

        SquaredSumIntensityNormalizer squaredSumIntensityNormalizer = new SquaredSumIntensityNormalizer();

        normalizedTestSpectra = new ISpectrum[testSpectra.size()];

        for (int i = 0; i < testSpectra.size(); i++) {
            normalizedTestSpectra[i] = new Spectrum(testSpectra.get(i), squaredSumIntensityNormalizer.normalizePeaks(testSpectra.get(i).getPeaks()));
        }
    }

//    @Test
//    public void testComparison() throws Exception {
//        final float PRECURSOR_TOLERANCE = 2;
//        final float FRAGMENT_TOLERANCE = 0.5F;
//
//        HypergeometricScore hypergeometricScore = new HypergeometricScore();
//        hypergeometricScore.setPeakMzTolerance(FRAGMENT_TOLERANCE);
//        ISimilarityChecker dotProduct = new FrankEtAlDotProduct(FRAGMENT_TOLERANCE, 15);
//        FisherExactTest fisherExactTest = new FisherExactTest();
//        fisherExactTest.setPeakMzTolerance(FRAGMENT_TOLERANCE);
//
//        IntensityRankCorrelation intensityRankCorrelation = new IntensityRankCorrelation(FRAGMENT_TOLERANCE);
//
//        BufferedWriter writer = new BufferedWriter(new FileWriter("/tmp/output" + FRAGMENT_TOLERANCE));
//
//        writer.write("Dot\thgt\tp\tint_rank\tprecursor_1\tprecursor_2\tid_1\tid_2\tidentical\tsame_seq_r1\tsame_seq_all_r\n");
//        for (int i = 0; i < testSpectra.size(); i++) {
//            for (int j = 0; j < testSpectra.size(); j++) {
//                if (j < i)
//                    continue;
//
//                ISpectrum spec1 = testSpectra.get(i);
//                ISpectrum spec2 = testSpectra.get(j);
//
//                if (Math.abs(spec1.getPrecursorMz() - spec2.getPrecursorMz()) > PRECURSOR_TOLERANCE)
//                    continue;
//
//                if ( i == 0 && j == 4)
//                    System.out.print("");
//
//                double hgt = hypergeometricScore.assessSimilarity(spec1, spec2);
//                double dot = dotProduct.assessSimilarity(spec1, spec2);
//                double p = fisherExactTest.assessSimilarity(spec1, spec2);
//                double intRank = intensityRankCorrelation.assessSimilarity(spec1, spec2);
//
//                String identical = (i == j) ? "TRUE" : "FALSE";
//                String spec1R1 = spec1.getProperty("splib_sequence");
//                String spec1R2 = spec1.getProperty("sequenceR2");
//                String spec2R1 = spec2.getProperty("splib_sequence");
//                String spec2R2 = spec2.getProperty("sequenceR2");
//                String sameSeq = spec1R1.equals(spec2R1) ? "TRUE" : "FALSE";
//                String sameSeqAll = (spec1R1.equals(spec2R1) || spec1R2.equals(spec2R1) || spec1R1.equals(spec2R2) || spec1R2.equals(spec2R2)) ? "TRUE" : "FALSE";
//
//                writer.write(dot + "\t" + hgt + "\t" + p + "\t" + intRank + "\t" + spec1.getPrecursorMz() + "\t"
//                        + spec2.getPrecursorMz() + "\t" + spec1.getId() + "\t"  + spec2.getId() + "\t" + identical + "\t" + sameSeq + "\t"  + sameSeqAll + "\n");
//            }
//        }
//
//        writer.close();
//    }

    @Test
    public void testPerformance() {
        /**
         * Results on my laptop (n = 50):
         * Dot Product: 1598 msec, 1142 msec, 1016 msec
         * Dot Product (getSharedPeaks): 1335 msec, 993 msec, 950 msec
         * Dot Product (getSharedPeaks2): 1269 msec, 1013 msec, 956 msec
         *
         * n = 5
         * HypergeometricScore: 11404 msec, 9308 msec, 9670 msec
         * Fisher Exact Test: Took: 4230 msec, 3570 msec, 3627 msec
         * Dot Product (new version): 432 msec, 109 msec, 108 msec   (300 - 500)
         *
         * -- changed CPU settings on my machine - slower now --
         * Combined version: 11-14 sec(!)
         * Fisher Exact: ~3500 msec
         * Hypergeometric: ~4300 msec
         * Dot: ~100-130 msec
         * Combined: 6500 msec
         *
         * -- CERN statistics --
         * rounds = 50, times = 4
         * Hypergeometric: ~7450 msec
         * Dot: ~820 msec
         * Combined: 3500 - 7000 msec (?)
         * Intensity Rank: 3000 msec
         */

        // generally, this test is not run
        if (true) return;

        int nRounds = 4;
        int nTimes = 4;

        ISimilarityChecker similarityChecker = new IntensityRankCorrelation(0.5F);
        long[] duration = new long[nTimes];

        for (int time = 0; time < nTimes; time++) {
            long startTime = System.currentTimeMillis();

            for (int round = 0; round < nRounds; round++) {
                for (int i = 0; i < testSpectra.size(); i++) {
                    for (ISpectrum aTestSpectra : testSpectra) {
                        similarityChecker.assessSimilarity(testSpectra.get(i), aTestSpectra);
                    }
                }
            }

            duration[time] = System.currentTimeMillis() - startTime;
        }

        for (long d : duration) {
            System.out.println("Took: " + d + " msec");
        }
    }
}

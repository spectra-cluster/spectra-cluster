package uk.ac.ebi.pride.spectracluster.similarity;

import junit.framework.TestCase;
import org.junit.Test;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.util.ClusteringTestUtilities;

import java.util.*;

/**
 * Created with IntelliJ IDEA.
 * User: jg
 * Date: 2/19/14
 * Time: 9:35 PM
 * To change this template use File | Settings | File Templates.
 */
public class ClosestPeaksTest {

    public static final String[] INTERESTING_IDS =
            {
                    "1248200",
                    "1249768",
                    "1249841",
                    "1250984",
//                    "1252317",
//                    "1252580",
//                    "1252923",
//                    "18834",
//                    "58604",
//                    "99609",
            };

    public static Set<String> INTERESTING_ID_SET = new HashSet<String>(Arrays.asList(INTERESTING_IDS));

    /**
     * really here to print out the peaks compared on a few interesting cases where the results are different
     */
    @Test
    public void showHandledPeaksForInterestingCases() {
        List<ISpectrum> spectra = ClusteringTestUtilities.readISpectraFromResource();

        ISpectrum[] spectrums = (ISpectrum[]) spectra.toArray();

        ISimilarityChecker checker = new FrankEtAlDotProductTester();
        ISimilarityChecker currentChecker = new FrankEtAlDotProduct(0.5, 15);

        //noinspection UnnecessaryLocalVariable,UnusedDeclaration,UnusedAssignment
        Set<String> interestingIds = new HashSet<String>();


        for (int i = 0; i < spectrums.length; i++) {
            ISpectrum psm1 = spectrums[i];
            String id1 = psm1.getId();
            if (!INTERESTING_ID_SET.contains(id1))
                continue; // not an interesting case

            for (int j = i + 1; j < spectrums.length; j++) {
                ISpectrum psm2 = spectrums[j];

                String id2 = psm2.getId();
                if (!INTERESTING_ID_SET.contains(id2))
                    continue; // not an interesting case

                //        System.out.println("Comparing " + id1 + " " + id2);

                StringBuilder usedPeaksTester = new StringBuilder();

                 //noinspection UnnecessaryLocalVariable,UnusedDeclaration,UnusedAssignment
                double dotOrg = checker.assessSimilarity(psm1, psm2);

                //         System.out.println("Peaks compared original Frank Et Al (when the code is written)");
                // print usage
                //        System.out.println(usedPeaksTester.toString());

                usedPeaksTester.setLength(0);  // clear debug output
                double dotNew = currentChecker.assessSimilarity(psm1, psm2);

                // print usage
                //          System.out.println("Peaks compared current Frank Et Al ");
                //          System.out.println(usedPeaksTester.toString());


            }

        }
    }

    /**
     * test that latest and refactored dot products give the same answer
     */
    // @Test  not sure this always works todo
    public void testDifferentDotProducts() {
        List<ISpectrum> spectra = ClusteringTestUtilities.readISpectraFromResource();

        ISpectrum[] spectrums = (ISpectrum[]) spectra.toArray();

        int total = 0;
        int different = 0;
        ISimilarityChecker checker = new FrankEtAlDotProduct(0.5, 15);
        ISimilarityChecker currentChecker = new FrankEtAlDotProductJohannes();

        Set<String> interestingIds = new HashSet<String>();


        for (int i = 0; i < spectrums.length; i++) {
            ISpectrum psm1 = spectrums[i];
            for (int j = i + 1; j < spectrums.length; j++) {
                ISpectrum psm2 = spectrums[j];
                double dotOrg = checker.assessSimilarity(psm1, psm2);
                double dotNew = currentChecker.assessSimilarity(psm1, psm2);

                if (Math.abs(dotOrg - dotNew) > 0.00001) {
                    different++;

                    StringBuilder usedPeaksTester = new StringBuilder();

                    // these are the really interesting cases
                    dotOrg = checker.assessSimilarity(psm1, psm2);

                    double noClosestPeak = dotNew;
                    dotNew = currentChecker.assessSimilarity(psm1, psm2);
                    String id2 = psm2.getId();
                    String id1 = psm1.getId();
                    interestingIds.add(id1);
                    interestingIds.add(id2);

                    //                 System.out.println(usedPeaksTester.toString());
                    //                 System.out.printf(id2 + ":" + id1 + " " + "Old: %8.3f Newx: %8.3f New: %8.3f\tDiff: %8.3f\n", dotOrg, noClosestPeak, dotNew, dotOrg - dotNew);
                }
                total++;

            }

        }

        List<String> sorted = new ArrayList<String>(interestingIds);
        Collections.sort(sorted);
        //      System.out.println("Interesting Ids");
        for (String s : sorted) {
            //         System.out.println(s);
        }


        TestCase.assertEquals(0, different);
    }

}

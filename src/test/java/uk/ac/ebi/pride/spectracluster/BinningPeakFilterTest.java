package uk.ac.ebi.pride.spectracluster;

import org.junit.*;
import uk.ac.ebi.pride.spectracluster.filter.*;
import uk.ac.ebi.pride.spectracluster.spectrum.*;
import uk.ac.ebi.pride.spectracluster.spectrum.Spectrum;
import uk.ac.ebi.pride.spectracluster.util.*;

import java.util.*;

/**
 * uk.ac.ebi.pride.spectracluster.BinningPeakFilterTest
 *
 * @author Steve Lewis
 * @date 29/05/2014
 */
public class BinningPeakFilterTest {

     @Test
    public void testBinningPeakFilterTest()
     {
         List<ISpectrum> spectra  = ClusteringTestUtilities.readISpectraFromResource();
         for (ISpectrum spectrum : spectra) {
             final List<IPeak> oldPeaks = spectrum.getPeaks();
             final List<IPeak> newPeaks = BinnedHighestNPeakFilter.DEFAULT.filter(oldPeaks);
             testFilteredPeaks(oldPeaks,newPeaks);
         }
     }

    private void testFilteredPeaks(List<IPeak> oldPeaks, List<IPeak> newPeaks) {
        for (int i = 0; i < 5000; i += 50) {
            testPeaksInBin(i,oldPeaks,newPeaks);

        }

    }

    public static final int BIN_SIZE =  BinnedHighestNPeakFilter.DEFAULT_BIN_SIZE;
    private void testPeaksInBin(int bin, List<IPeak> oldPeaks, List<IPeak> newPeaks) {
        int oldCount = 0;
        int newCount = 0;
        double endBin = BIN_SIZE + bin;
        for (IPeak peak : oldPeaks) {
             double mz =  peak.getMz();
             if(mz < bin)
                 continue;
            if(mz > endBin)
                break;
            oldCount++;
        }
        for (IPeak peak : newPeaks) {
             double mz =  peak.getMz();
             if(mz < bin)
                 continue;
            if(mz > endBin)
                break;
            newCount++;
        }

        oldCount = Math.min(oldCount,BinnedHighestNPeakFilter.DEFAULT_MAX_PEAKS_PER_BIN);
        if(newCount < oldCount) {
            if(newCount < oldCount - 1) {
                if(newCount < oldCount - 2) {
                    if(newCount < oldCount - 3) {
                        Assert.assertTrue(newCount >= oldCount);
                      }
                  }
              }
        }

    }
}

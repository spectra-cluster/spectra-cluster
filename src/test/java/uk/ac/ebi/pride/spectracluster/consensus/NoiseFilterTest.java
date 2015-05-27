package uk.ac.ebi.pride.spectracluster.consensus;

import junit.framework.Assert;
import org.junit.Before;
import org.junit.Test;
import uk.ac.ebi.pride.spectracluster.spectrum.IPeak;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.util.ClusteringTestUtilities;
import uk.ac.ebi.pride.spectracluster.util.MZIntensityUtilities;
import uk.ac.ebi.pride.spectracluster.util.comparator.PeakIntensityComparator;
import uk.ac.ebi.pride.spectracluster.util.comparator.PeakMzComparator;
import uk.ac.ebi.pride.spectracluster.util.function.IFunction;
import uk.ac.ebi.pride.spectracluster.util.function.peak.BinnedHighestNPeakFunction;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Created by jg on 16.12.14.
 */
public class NoiseFilterTest {
    private List<ISpectrum> testSpectra;
    public final int PEAKS_TO_KEEP = 5;

    @Before
    public void setUp() {
        testSpectra = ClusteringTestUtilities.readISpectraFromResource();
    }

    @Test
    public void testBinnedHighestNPeakFilter() {
        IFunction<List<IPeak>, List<IPeak>> testFilter = new BinnedHighestNPeakFunction(PEAKS_TO_KEEP, (int) ConsensusSpectrum.NOISE_FILTER_INCREMENT, 0);

        for (ISpectrum spectrum : testSpectra) {
            List<IPeak> originalPeaks = spectrum.getPeaks();

            List<IPeak> classicalFilterPeaks = classicalFilterNoise(originalPeaks, PEAKS_TO_KEEP);
            List<IPeak> testFilterPeaks = testFilter.apply(originalPeaks);

            Assert.assertEquals(classicalFilterPeaks.size(), testFilterPeaks.size());

            Collections.sort(testFilterPeaks);

            for (int i = 0; i < classicalFilterPeaks.size(); i++) {
                IPeak classicalFilterPeak = classicalFilterPeaks.get(i);
                IPeak testFilterPeak = testFilterPeaks.get(i);

                Assert.assertEquals(classicalFilterPeak, testFilterPeak);
            }
        }
    }

    /**
     * This function is the original peak filtering function extracted from the
     * ConsensusSpectrum class.
     * @param inp
     * @param peaksInBinToKeep
     * @return
     */
    protected static List<IPeak> classicalFilterNoise(List<IPeak> inp, int peaksInBinToKeep) {
        List<IPeak> filteredSpectrum = new ArrayList<IPeak>();

        int lowerBound = 0;
        // process the peaks using a sliding window of 100 m/z
        for (double startMz = 0, endMz = ConsensusSpectrum.NOISE_FILTER_INCREMENT; endMz <= MZIntensityUtilities.HIGHEST_USABLE_MZ; endMz += ConsensusSpectrum.NOISE_FILTER_INCREMENT, startMz += ConsensusSpectrum.NOISE_FILTER_INCREMENT) {
            List<IPeak> peakBuffer = new ArrayList<IPeak>();

            // set the lower bound
            for (int i = lowerBound; i < inp.size(); i++) {
                if (inp.get(i).getMz() >= startMz) {
                    lowerBound = i;
                    break;
                }
            }

            if (inp.get(lowerBound).getMz() < startMz)
                continue;

            for (int i = lowerBound; i < inp.size(); i++) {
                if (inp.get(i).getMz() <= endMz) {
                    peakBuffer.add(inp.get(i));
                } else {
                    lowerBound = i;
                    break;
                }
            }

            if (peakBuffer.size() < 1)
                continue;

            Collections.sort(peakBuffer, PeakIntensityComparator.INSTANCE);

            List<IPeak> peaksInBin = new ArrayList<IPeak>(peaksInBinToKeep);

            for (int i = 0; i < peaksInBinToKeep && i < peakBuffer.size(); i++)
                peaksInBin.add(peakBuffer.get(i));

            Collections.sort(peaksInBin, new PeakMzComparator());
            filteredSpectrum.addAll(peaksInBin);
        }

        return filteredSpectrum;
    }
}

package uk.ac.ebi.pride.spectracluster.consensus;

import junit.framework.Assert;
import org.junit.Test;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.util.ClusteringTestUtilities;
import uk.ac.ebi.pride.spectracluster.util.Defaults;
import uk.ac.ebi.pride.spectracluster.util.JMZTabUtilities;
import uk.ac.ebi.pride.tools.jmzreader.model.Spectrum;
import uk.ac.ebi.pride.tools.pride_spectra_clustering.util.Peak;

import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * @author Rui Wang
 * @version $Id$
 */
public class MergeIdenticalPeaksComparisonTests {


    @Test
    public void testSpectrumBuilderComparison() throws Exception {

        //noinspection UnusedDeclaration
        IConsensusSpectrumBuilder frankEtAlConsensusSpectrumBuilder = Defaults.getDefaultConsensusSpectrumBuilder();
        //noinspection UnusedDeclaration
        uk.ac.ebi.pride.tools.pride_spectra_clustering.consensus_spectrum_builder.impl.FrankEtAlConsensusSpectrumBuilder        //noinspection UnusedDeclaration
                originalFrankEtAlConsensusSpectrumBuilder = new uk.ac.ebi.pride.tools.pride_spectra_clustering.consensus_spectrum_builder.impl.FrankEtAlConsensusSpectrumBuilder();

        List<ISpectrum> spectra = ClusteringTestUtilities.readISpectraFromResource();
        List<Spectrum> originalSpectra = JMZTabUtilities.readSpectrumsFromResource();

        Assert.assertEquals(spectra.size(), originalSpectra.size());

        // todo rewrite
//        for (int i = 0; i < spectra.size() - 1; i++) {
//            ISpectrum spectrum = spectra.get(i);
//            ISpectrum nextSpectrum = spectra.get(i + 1);
//
//            List<IPeak> peaks = frankEtAlConsensusSpectrumBuilder.addAllPeaks(Arrays.asList(spectrum, nextSpectrum));
//            List<Peak> originalPeaks = Adapters.fromIPeaks(peaks);
//            Map<Double, Peak> originalPeakMap = createPeakMap(originalPeaks);
//            List<IPeak> mergedPeaks = frankEtAlConsensusSpectrumBuilder.mergeIdenticalPeaks(peaks);
//            List<Peak> originalMergedPeaks = originalFrankEtAlConsensusSpectrumBuilder.mergeIdenticalPeaks(originalPeakMap);
//            List<IPeak> convertedOriginalMergedPeaks = Adapters.fromPeaks(originalMergedPeaks);
//
//            Collections.sort(mergedPeaks);
//            Collections.sort(convertedOriginalMergedPeaks);
//
//            Assert.assertEquals(mergedPeaks.size(), convertedOriginalMergedPeaks.size());
//
//            for (int j = 0; j < mergedPeaks.size(); j++) {
//                IPeak mergedPeak = mergedPeaks.get(j);
//                IPeak originalMergedPeak = convertedOriginalMergedPeaks.get(j);
//                boolean equivalent = mergedPeak.equivalent(originalMergedPeak);
//                if (!equivalent) {
//                    Assert.assertTrue(equivalent);
//                }
//            }
//
//        }
//
    }


    @SuppressWarnings("UnusedDeclaration")
    private Map<Double, Peak> createPeakMap(List<Peak> originalPeaks) {
        Map<Double, Peak> originalPeakMap = new LinkedHashMap<Double, Peak>();

        for (Peak originalPeak : originalPeaks) {
            originalPeakMap.put(originalPeak.getMz(), originalPeak);
        }


        return originalPeakMap;
    }

}

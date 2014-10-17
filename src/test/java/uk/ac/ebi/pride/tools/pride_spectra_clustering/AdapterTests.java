package uk.ac.ebi.pride.tools.pride_spectra_clustering;

import org.junit.*;
import uk.ac.ebi.pride.spectracluster.spectrum.IPeak;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.util.*;
import uk.ac.ebi.pride.tools.jmzreader.model.Spectrum;
import uk.ac.ebi.pride.tools.pride_spectra_clustering.impl.Adapters;
import uk.ac.ebi.pride.tools.pride_spectra_clustering.util.ClusteringSpectrum;
import uk.ac.ebi.pride.tools.pride_spectra_clustering.util.Peak;

import java.util.List;

/**
 * uk.ac.ebi.pride.tools.pride_spectra_clustering.AdapterTests
 * User: Steve
 * Date: 6/19/13
 */
public class AdapterTests {


    @Test
    public void testSpectra() throws Exception {
        List<ISpectrum> sprctra = ClusteringTestUtilities.readConsensusSpectralItems();
        final List<Spectrum> spectrums = JMZTabUtilities.readSpectrumsFromResource();

        Assert.assertEquals(sprctra.size(), spectrums.size());
        for (int i = 0; i < sprctra.size(); i++) {
            ISpectrum readISpectrum = sprctra.get(i);
            Spectrum readSpectrum = spectrums.get(i);

            ClusteringSpectrum readClusteringSpectrum = new ClusteringSpectrum(readSpectrum);

            ISpectrum adapted = Adapters.fromSpectrum(readClusteringSpectrum);
            Assert.assertTrue(readISpectrum.equivalent(adapted));

        }

    }

    @Test
    public void testPeaks() throws Exception {
        List<ISpectrum> sprctra = ClusteringTestUtilities.readConsensusSpectralItems();
        final List<Spectrum> spectrums = JMZTabUtilities.readSpectrumsFromResource();

        Assert.assertEquals(sprctra.size(), spectrums.size());
        for (int i = 0; i < sprctra.size(); i++) {
            ISpectrum readISpectrum = sprctra.get(i);
            final List<IPeak> peaks = readISpectrum.getPeaks();

            Spectrum readSpectrum = spectrums.get(i);
            ClusteringSpectrum readClusteringSpectrum = new ClusteringSpectrum(readSpectrum);
            final List<Peak> peaklist = readClusteringSpectrum.getPeaklist();
            List<IPeak> concertedpeaks = Adapters.fromPeaks(peaklist);
            Assert.assertEquals(peaks.size(), concertedpeaks.size());
            for (int j = 0; j < peaklist.size(); j++) {
                Assert.assertTrue(peaks.get(j).equivalent(concertedpeaks.get(j)));
            }

            ISpectrum adapted = Adapters.fromSpectrum(readClusteringSpectrum);
            boolean equivalent = readISpectrum.equivalent(adapted);
            if(!equivalent)
                equivalent = readISpectrum.equivalent(adapted); // break here
            Assert.assertTrue(equivalent);

        }

    }
}

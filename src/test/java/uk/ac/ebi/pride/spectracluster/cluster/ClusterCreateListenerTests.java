package uk.ac.ebi.pride.spectracluster.cluster;


import org.junit.*;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.util.ClusteringTestUtilities;
import uk.ac.ebi.pride.spectracluster.io.ParserUtilities;
import uk.ac.ebi.pride.spectracluster.util.SpectrumCreateListener;

import java.io.File;
import java.io.FileReader;
import java.io.LineNumberReader;
import java.util.ArrayList;
import java.util.List;

/**
 * uk.ac.ebi.pride.spectracluster.cluster.ClusterCreateListenerTests
 *
 * @author Steve Lewis
 * @date 29/10/13
 */
public class ClusterCreateListenerTests {

    @Test
    public void testSpectrumListeners() throws Exception {
        List<ISpectrum> originalSpectra = ClusteringTestUtilities.readConsensusSpectralItems();
        File spectrumFile = ClusteringTestUtilities.getSpectrumFile();
        LineNumberReader rdr = new LineNumberReader(new FileReader(spectrumFile));

        TestClusterListener listener = new TestClusterListener();
        ParserUtilities.readAndProcessSpectra(rdr, listener);

        List<ISpectrum> readSpectra = listener.getSpectra();

        Assert.assertEquals(readSpectra.size(), originalSpectra.size());

        for (int i = 0; i < readSpectra.size(); i++) {
            ISpectrum os = originalSpectra.get(i);
            ISpectrum ls = readSpectra.get(i);
            Assert.assertTrue(os.equivalent(ls));

        }

    }

    public static class TestClusterListener implements SpectrumCreateListener {

        private final List<ISpectrum> m_Spectra = new ArrayList<ISpectrum>();


        public List<ISpectrum> getSpectra() {
            return new ArrayList<ISpectrum>(m_Spectra);
        }

        /**
         * initialize reading - if reading happens once - say  from
         * one file all this may happen in the constructor
         */
        @Override
        public void onSpectrumStarted() {
            m_Spectra.clear();

        }

        /**
         * do something when a Spectrum is created or read
         *
         * @param Spectrum
         */
        @Override
        public void onSpectrumCreate(ISpectrum spec) {
            m_Spectra.add(spec);

        }

        /**
         * do something when a Spectrum when the last Spectrum is read -
         * this may be after a file read is finished
         */
        @Override
        public void onSpectrumCreateFinished() {
            // nothing to do

        }

    }
}

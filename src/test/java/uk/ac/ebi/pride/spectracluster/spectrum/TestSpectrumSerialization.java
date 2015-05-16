package uk.ac.ebi.pride.spectracluster.spectrum;

import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;
import uk.ac.ebi.pride.spectracluster.io.ParserUtilities;

import java.io.*;
import java.util.ArrayList;
import java.util.List;

/**
 * Created by jg on 13.05.15.
 */
public class TestSpectrumSerialization {
    private List<Spectrum> testSpectra;

    @Before
    public void setUp() throws Exception {
        File testFile = new File(TestSpectrumSerialization.class.getClassLoader().getResource("spectra_400.0_4.0.mgf").toURI());
        testSpectra = new ArrayList<Spectrum>();
        ISpectrum[] readSpectra = ParserUtilities.readMGFScans(testFile);

        for (ISpectrum s : readSpectra)
            testSpectra.add((Spectrum) s);
    }

    @Test
    public void testSerialization() throws Exception {
        for (int i = 0; i < 10; i++) {
            Spectrum s = testSpectra.get(i);

            ByteArrayOutputStream outputStream = new ByteArrayOutputStream();
            ObjectOutputStream objectOutputStream = new ObjectOutputStream(outputStream);

            objectOutputStream.writeObject(s);

            objectOutputStream.close();
            outputStream.close();

            ByteArrayInputStream byteInputStream = new ByteArrayInputStream(outputStream.toByteArray());
            ObjectInputStream inputStream = new ObjectInputStream(byteInputStream);

            Spectrum sLoaded = (Spectrum) inputStream.readObject();

            Assert.assertEquals(s.getId(), sLoaded.getId());
            Assert.assertEquals(s.getPrecursorCharge(), sLoaded.getPrecursorCharge());
            Assert.assertEquals(s.getPrecursorMz(), sLoaded.getPrecursorMz(), 0);
            Assert.assertEquals(s.getSumSquareIntensity(), sLoaded.getSumSquareIntensity(), 0);
            Assert.assertEquals(s.getPeaks().size(), sLoaded.getPeaks().size());

            for (int j = 0; j < s.getPeaks().size(); j++) {
                Assert.assertEquals(s.getPeaks().get(j).getMz(), sLoaded.getPeaks().get(j).getMz(), 0);
                Assert.assertEquals(s.getPeaks().get(j).getIntensity(), sLoaded.getPeaks().get(j).getIntensity(), 0);
                Assert.assertEquals(s.getPeaks().get(j).getCount(), sLoaded.getPeaks().get(j).getCount(), 0);
            }
        }
    }

    @Test
    public void testAsMajorPeaks() {
        ISpectrum spectrum = testSpectra.get(0);

        int[] majorPeaks = spectrum.asMajorPeakMZs(10);

        Assert.assertEquals(10, majorPeaks.length);
        Assert.assertEquals(199, majorPeaks[0]);
        Assert.assertEquals(227, majorPeaks[1]);
        Assert.assertEquals(274, majorPeaks[2]);
        Assert.assertEquals(345, majorPeaks[3]);
        Assert.assertEquals(454, majorPeaks[4]);
        Assert.assertEquals(514, majorPeaks[5]);
        Assert.assertEquals(525, majorPeaks[6]);
        Assert.assertEquals(548, majorPeaks[7]);
        Assert.assertEquals(572, majorPeaks[8]);
        Assert.assertEquals(685, majorPeaks[9]);
    }
}

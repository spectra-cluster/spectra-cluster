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
}

package uk.ac.ebi.pride.spectracluster;

import org.junit.*;
import uk.ac.ebi.pride.spectracluster.spectrum.IPeak;
import uk.ac.ebi.pride.spectracluster.spectrum.Peak;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

/**
 * uk.ac.ebi.pride.spectracluster.SpectrumPeakTests
 *
 * @author Steve Lewis
 */
public class SpectrumPeakTests {

    public static final Random RND = new Random();
    public static final int MAX_PEAKS = 100;


    public static IPeak[] buildPeaks() {
        List<IPeak> holder = new ArrayList<>();
        for (int i = 0; i < MAX_PEAKS; i++) {
            float mz = 1000 - i;
            float intensity = 100 * RND.nextFloat();
            holder.add(new Peak(mz, intensity));

        }

        IPeak[] ret = new IPeak[holder.size()];
        holder.toArray(ret);
        return ret;
    }

    @Test
    public void sortTest() {
        IPeak[] pks = buildPeaks();
        Assert.assertEquals(MAX_PEAKS, pks.length);
        double last = Double.MAX_VALUE;
        // we built them in inverted order
        for (IPeak pk : pks) {
            Assert.assertTrue(pk.getMz() < last);
            last = pk.getMz();
        }
        Arrays.sort(pks);
        last = Double.MIN_VALUE;
        // they sort in ascendint order
        for (IPeak pk : pks) {
            Assert.assertTrue(pk.getMz() > last);
            last = pk.getMz();
        }
    }

}

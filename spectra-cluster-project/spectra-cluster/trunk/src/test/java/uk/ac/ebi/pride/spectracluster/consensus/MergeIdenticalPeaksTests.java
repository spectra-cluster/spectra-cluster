package uk.ac.ebi.pride.spectracluster.consensus;

import org.junit.*;
import uk.ac.ebi.pride.spectracluster.spectrum.IPeak;
import uk.ac.ebi.pride.spectracluster.spectrum.Peak;
import uk.ac.ebi.pride.spectracluster.util.*;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Rui Wang
 * @version $Id$
 */
public class MergeIdenticalPeaksTests {

    @SuppressWarnings({"FieldCanBeLocal", "MismatchedQueryAndUpdateOfCollection"})
    private List<IPeak> peaks;
    @SuppressWarnings({"FieldCanBeLocal", "UnusedDeclaration"})
    private IConsensusSpectrumBuilder consensusSpectrumBuilder;

    @Before
    public void setUp() throws Exception {
        peaks = new ArrayList<IPeak>();

        peaks.add(new Peak(123.45F, 10, 1));
        peaks.add(new Peak(123.47F, 20, 1));
        peaks.add(new Peak(123.48F, 30, 1));

        peaks.add(new Peak(223.0F, 10, 1));
        peaks.add(new Peak(223.1F, 20, 1));
        peaks.add(new Peak(223.4F, 30, 1));

        consensusSpectrumBuilder = Defaults.getDefaultConsensusSpectrumBuilder();
    }

    @Test
    public void testIdenticalPeaksMerged() throws Exception {
//        List<IPeak> mergedPeaks = consensusSpectrumBuilder.mergeIdenticalPeaksInternal(peaks);
//
//        Assert.assertEquals(4, mergedPeaks.size());
//        Assert.assertEquals(60.0F, mergedPeaks.get(0).getIntensity());
//        Assert.assertEquals(60.0F, mergedPeaks.get(1).getIntensity());
    }
}

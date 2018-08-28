package uk.ac.ebi.pride.spectracluster.cluster;

import junit.framework.Assert;
import org.junit.Before;
import org.junit.Test;
import uk.ac.ebi.pride.spectracluster.quality.OriginalSignalToNoiseChecker;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.util.ClusterUtilities;
import uk.ac.ebi.pride.spectracluster.util.ClusteringTestUtilities;
import uk.ac.ebi.pride.spectracluster.util.comparator.DefaultClusterComparator;

import java.util.List;

/**
 * @author Rui Wang
 * @version $Id$
 */
public class ClusterComparatorTests {

    private List<ISpectrum> peptideSpectrumMatches;
    private DefaultClusterComparator defaultClusterComparator;
    private OriginalClusterComparator originalClusterComparator;

    @Before
    public void setUp() throws Exception {
        defaultClusterComparator = new DefaultClusterComparator();
        originalClusterComparator = new OriginalClusterComparator(new OriginalSignalToNoiseChecker());

        peptideSpectrumMatches = ClusteringTestUtilities.readConsensusSpectralItems();
    }

    @Test
    public void testClusterComparator() throws Exception {
        final ICluster cluster1 = ClusterUtilities.asCluster(peptideSpectrumMatches.get(0));
        int result = defaultClusterComparator.compare(cluster1, ClusterUtilities.asCluster(peptideSpectrumMatches.get(1)));
        int originalResult = originalClusterComparator.compare(cluster1, ClusterUtilities.asCluster(peptideSpectrumMatches.get(1)));
        Assert.assertEquals(result, originalResult);
    }
}

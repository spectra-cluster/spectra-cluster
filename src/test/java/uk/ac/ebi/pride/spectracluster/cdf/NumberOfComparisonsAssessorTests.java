package uk.ac.ebi.pride.spectracluster.cdf;

import junit.framework.Assert;
import org.junit.Test;
import uk.ac.ebi.pride.spectracluster.cluster.ICluster;
import uk.ac.ebi.pride.spectracluster.cluster.SpectralCluster;
import uk.ac.ebi.pride.spectracluster.consensus.ConsensusSpectrum;
import uk.ac.ebi.pride.spectracluster.spectrum.IPeak;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Created by jg on 13.10.17.
 */
public class NumberOfComparisonsAssessorTests {
    @Test
    public void testMinNumberComparisonAssessor() {
        INumberOfComparisonAssessor assessor = new MinNumberComparisonsAssessor(5000);

        Assert.assertEquals(5000, assessor.getNumberOfComparisons(null, 5));
        Assert.assertEquals(7000, assessor.getNumberOfComparisons(null, 7000));
    }

    @Test
    public void testSpectraPerBinNumberComparisonsAssessor() {
        SpectraPerBinNumberComparisonAssessor assessor = new SpectraPerBinNumberComparisonAssessor(2F);

        float[] existingPrecursors = {1F, 2F, 3F, 4F, 5F, 6F, 10F, 11F, 15F, 18F};
        for (float prec : existingPrecursors) {
            assessor.countSpectrum(prec);

        }

        List<IPeak> peaklist = new ArrayList<IPeak>();

        ICluster cluster = new SpectralCluster("test", new ConsensusSpectrum("est", 1, 2,
                2, 2, peaklist, 0.0F));

        Assert.assertEquals(2, assessor.getNumberOfComparisons(cluster, 4));

        cluster = new SpectralCluster("test", new ConsensusSpectrum("est", 1, 18,
                2, 2, peaklist, 0.0F));
        Assert.assertEquals(1, assessor.getNumberOfComparisons(cluster, 1));
    }
}

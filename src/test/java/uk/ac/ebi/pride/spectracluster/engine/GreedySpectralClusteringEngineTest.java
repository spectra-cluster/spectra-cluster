package uk.ac.ebi.pride.spectracluster.engine;

import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;
import uk.ac.ebi.pride.spectracluster.cluster.ICluster;
import uk.ac.ebi.pride.spectracluster.io.ParserUtilities;
import uk.ac.ebi.pride.spectracluster.similarity.CombinedFisherIntensityTest;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.util.ClusterUtilities;
import uk.ac.ebi.pride.spectracluster.util.Defaults;
import uk.ac.ebi.pride.spectracluster.util.function.peak.FractionTICPeakFunction;

import java.io.File;
import java.util.*;

/**
 * Created by jg on 11.05.15.
 */
public class GreedySpectralClusteringEngineTest {
    private List<ISpectrum> testSpectra;

    @Before
    public void setUp() throws Exception {
        File testFile = new File(GreedySpectralClusteringEngineTest.class.getClassLoader().getResource("spectra_400.0_4.0.mgf").toURI());
        testSpectra = new ArrayList<>();
        ISpectrum[] readSpectra = ParserUtilities.readMGFScans(testFile);

        testSpectra.addAll(Arrays.asList(readSpectra));

        testSpectra.sort(new SpectrumMzComparator());
    }

    @Test
    public void testClustering() throws Exception {
        List<ICluster> cluster = new ArrayList<>();
        List<ICluster> secondClustes = new ArrayList<>();

        GreedyIncrementalClusteringEngine engine = new GreedyIncrementalClusteringEngine(new CombinedFisherIntensityTest(0.5F), Defaults.getDefaultSpectrumComparator(), 4F, 0.6, new FractionTICPeakFunction(0.5F, 20));
        GreedyIncrementalClusteringEngine secondEngine = new GreedyIncrementalClusteringEngine(new CombinedFisherIntensityTest(0.5F), Defaults.getDefaultSpectrumComparator(), 4F, 0.95, new FractionTICPeakFunction(0.5F, 20));

        for (ISpectrum s : testSpectra) {
            ICluster spectrumCluster = ClusterUtilities.asCluster(s);

            // use engine
            List<ICluster> removedClusters = engine.addClusterIncremental(spectrumCluster);
            cluster.addAll(removedClusters);

            // use second engine
            List<ICluster> secondRemovedClusters = secondEngine.addClusterIncremental(spectrumCluster);
            secondClustes.addAll(secondRemovedClusters);
        }

        Assert.assertEquals(0, cluster.size());
        Assert.assertEquals(1, engine.getClusters().size());

        Assert.assertEquals(0, secondClustes.size());
        Assert.assertEquals(12, secondEngine.getClusters().size());
    }

    public class SpectrumMzComparator implements Comparator<ISpectrum> {
        @Override
        public int compare(ISpectrum o1, ISpectrum o2) {
            return Float.compare(o1.getPrecursorMz(), o2.getPrecursorMz());
        }
    }
}

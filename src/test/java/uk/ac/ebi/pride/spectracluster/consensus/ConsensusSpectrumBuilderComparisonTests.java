package uk.ac.ebi.pride.spectracluster.consensus;

import junit.framework.Assert;
import org.junit.Test;
import uk.ac.ebi.pride.spectracluster.filter.IPeakFilter;
import uk.ac.ebi.pride.spectracluster.filter.MaximialPeakFilter;
import uk.ac.ebi.pride.spectracluster.similarity.AllPeaksDotProduct;
import uk.ac.ebi.pride.spectracluster.spectrum.IPeak;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.spectrum.Spectrum;
import uk.ac.ebi.pride.spectracluster.util.ClusteringTestUtilities;
import uk.ac.ebi.pride.spectracluster.util.Defaults;

import java.util.ArrayList;
import java.util.List;

/**
 * Unit tests on ClusteringEngine to compare the clusters produced by AllPeaksDotProduct and
 * FrankAtElDotProduct.
 *
 *
 *
 * @author Rui Wang
 * @version $Id$
 */
public class ConsensusSpectrumBuilderComparisonTests {
    public static final String SMALL_CLUSTER_FILE = "uk/ac/ebi/pride/spectracluster/engine/cluster_spec_39546.mgf";

    @Test
    public void testSpectraFiltering() throws Exception {
        // filter the spectra first
        IPeakFilter filter = new MaximialPeakFilter(MaximialPeakFilter.DEFAULT_MAX_PEAKS);

        // load and filter the spectra
        List<ISpectrum> spectra= ClusteringTestUtilities.readISpectraFromResource(SMALL_CLUSTER_FILE);
        List<ISpectrum> filteredSpectra = new ArrayList<ISpectrum>(spectra.size());

        for (ISpectrum s : spectra) {
            List<IPeak> peaks = s.getPeaks();
            List<IPeak> filteredPeaks = filter.filter(peaks);

            Spectrum filteredSpectrum = new Spectrum(s, filteredPeaks);

            filteredSpectra.add(filteredSpectrum);
        }

        // create the consensus spectra
        final IConsensusSpectrumBuilder builder = Defaults.getDefaultConsensusSpectrumBuilder();

        builder.addSpectra(spectra.toArray(new ISpectrum[spectra.size()]));
        ISpectrum unfilteredConsensusSpectrum = builder.getConsensusSpectrum();

        builder.clear();
        builder.addSpectra(filteredSpectra.toArray(new ISpectrum[filteredSpectra.size()]));
        ISpectrum filteredConsensusSpectrum = builder.getConsensusSpectrum();

        // compare the spectra
        final AllPeaksDotProduct allPeaksDotProduct = new AllPeaksDotProduct(Defaults.getSimilarityMZRange());

        final double similarity = allPeaksDotProduct.assessSimilarity(unfilteredConsensusSpectrum, filteredConsensusSpectrum);
        Assert.assertTrue("Consensus spectra must be similar (" + similarity + ")", similarity >= 0.92);
    }
}

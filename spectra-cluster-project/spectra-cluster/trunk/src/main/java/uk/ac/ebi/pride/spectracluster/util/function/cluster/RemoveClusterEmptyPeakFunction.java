package uk.ac.ebi.pride.spectracluster.util.function.cluster;

import uk.ac.ebi.pride.spectracluster.cluster.ICluster;
import uk.ac.ebi.pride.spectracluster.cluster.SpectralCluster;
import uk.ac.ebi.pride.spectracluster.consensus.IConsensusSpectrumBuilder;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.util.function.IFunction;
import uk.ac.ebi.pride.spectracluster.util.function.spectrum.RemoveSpectrumEmptyPeakFunction;

/**
 * Remove all empty peaks from a cluster, including both the clustered spectra and consensus spectrum
 *
 * @author Rui Wang
 * @version $Id$
 */
public class RemoveClusterEmptyPeakFunction implements IFunction<ICluster, ICluster> {
    private final RemoveSpectrumEmptyPeakFunction removeSpectrumEmptyPeakFunction;

    public RemoveClusterEmptyPeakFunction() {
        removeSpectrumEmptyPeakFunction = new RemoveSpectrumEmptyPeakFunction();
    }

    @Override
    public ICluster apply(ICluster cluster) {
        IConsensusSpectrumBuilder consensusSpectrumBuilder = cluster.getConsensusSpectrumBuilder();
        String clusterId = cluster.getId();

        SpectralCluster processedCluster = new SpectralCluster(clusterId, consensusSpectrumBuilder);
        for (ISpectrum spectrum : cluster.getClusteredSpectra()) {
            ISpectrum vettedSpectrum = removeSpectrumEmptyPeakFunction.apply(spectrum);
            processedCluster.addSpectra(vettedSpectrum);
        }

        return processedCluster;
    }
}

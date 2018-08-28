package uk.ac.ebi.pride.tools.pride_spectra_clustering.impl;

import uk.ac.ebi.pride.spectracluster.cluster.ICluster;
import uk.ac.ebi.pride.spectracluster.cluster.SpectralCluster;
import uk.ac.ebi.pride.spectracluster.spectrum.IPeak;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.spectrum.KnownProperties;
import uk.ac.ebi.pride.spectracluster.spectrum.Spectrum;
import uk.ac.ebi.pride.spectracluster.util.Defaults;
import uk.ac.ebi.pride.tools.jmzreader.model.impl.ParamGroup;
import uk.ac.ebi.pride.tools.pride_spectra_clustering.consensus_spectrum_builder.ConsensusSpectrumBuilder;
import uk.ac.ebi.pride.tools.pride_spectra_clustering.consensus_spectrum_builder.impl.FrankEtAlConsensusSpectrumBuilder;
import uk.ac.ebi.pride.tools.pride_spectra_clustering.normalizer.IntensityNormalizer;
import uk.ac.ebi.pride.tools.pride_spectra_clustering.normalizer.impl.TotalIntensityNormalizer;
import uk.ac.ebi.pride.tools.pride_spectra_clustering.util.ClusteringSpectrum;
import uk.ac.ebi.pride.tools.pride_spectra_clustering.util.Peak;
import uk.ac.ebi.pride.tools.pride_spectra_clustering.util.SpectraCluster;

import java.util.Collections;
import java.util.List;
import java.util.Objects;
import java.util.stream.Collectors;

/**
 * uk.ac.ebi.pride.tools.pride_spectra_clustering.impl.Adapters
 * Static methods to convert to and from old code structures to new
 * User: Steve
 * Date: 6/18/13
 */
public class Adapters {

    public static final ConsensusSpectrumBuilder consensusBuilder = new FrankEtAlConsensusSpectrumBuilder();
    /**
     * The intensity normalizer used
     */
    public static final IntensityNormalizer normalizer = new TotalIntensityNormalizer();


    /**
     * convert old  SpectraCluster to new
     *
     * @param inp !null cluster
     * @return !null equivalent cluster
     */
    public static ICluster fromSpectraCluster(SpectraCluster inp) {
        String id = "";
        ICluster ret = new SpectralCluster(id,Defaults.getDefaultConsensusSpectrumBuilder());
        final List<ClusteringSpectrum> spectra = inp.getSpectra();

        for (ClusteringSpectrum sc : spectra) {
            ISpectrum spc = fromSpectrum(sc);
            ret.addSpectra(spc);
        }
        return ret;
    }

    /**
     * convert old  Spectrum to new
     *
     * @param inp !null Spectrum
     * @return !null equivalent Spectrum
     */
    public static ISpectrum fromSpectrum(ClusteringSpectrum inp) {
        final Integer precursorCharge = inp.getPrecursorCharge();
        final double precursorMZ = inp.getPrecursorMZ();
        final String id = inp.getId();
        final ParamGroup additional = inp.getAdditional();
        String peptide = inp.getPeptide();
        final List<Peak> peaklist = inp.getPeaklist();
        List<IPeak> newPeaks = peaklist.stream().map(Adapters::fromPeak).sorted().collect(Collectors.toList());

        ISpectrum ret = new Spectrum(id,
                 (int) (precursorCharge + 0.5),
                (float) precursorMZ,
                 Defaults.getDefaultQualityScorer(),
                newPeaks);

        ret.setProperty(KnownProperties.ANNOTATION_KEY,peptide);
        ret.setProperty(KnownProperties.IDENTIFIED_PEPTIDE_KEY,peptide);
         return ret;
    }

    /**
     * convert old peak to new  peak
     *
     * @param pk !null peak
     * @return !null equivalent peak
     */
    public static IPeak fromPeak(Peak pk) {
        return new uk.ac.ebi.pride.spectracluster.spectrum.Peak((float) pk.getMz(), (float) pk.getIntensity(), pk.getCount());
    }


    /**
     * convert old peak to new  peak
     *
     * @return !null equivalent peak list
     */
    public static List<IPeak> fromPeaks(List<Peak> pks) {
        List<IPeak> ret = pks.stream().filter(Objects::nonNull).map(Adapters::fromPeak).sorted().collect(Collectors.toList());
        return ret;
    }

    /**
     * convert old peak to new  peak
     *
     * @return !null equivalent peak list
     */
    public static List<Peak> fromIPeaks(List<IPeak> pks) {
        Collections.sort(pks);
        List<Peak> ret = pks.stream().map(Adapters::fromPeak).collect(Collectors.toList());
        return ret;
    }


    /**
     * convert old  SpectraCluster to new
     *
     * @param inp !null cluster
     * @return !null equivalent cluster
     */
    public static SpectraCluster fromSpectraCluster(ICluster inp) {
        String id = "";
        SpectraCluster ret = new SpectraCluster(consensusBuilder, normalizer);
        final List<ISpectrum> spectra = inp.getClusteredSpectra();

        for (ISpectrum sc : spectra) {
            ClusteringSpectrum spc = fromSpectrum(sc);
            ret.addSpectrum(spc);
        }
        return ret;
    }


    /**
     * convert new  Spectrum to old
     *
     * @param inp !null Spectrum
     * @return !null equivalent Spectrum
     */
    public static ClusteringSpectrum fromSpectrum(ISpectrum inp) {
        final double precursorCharge = inp.getPrecursorCharge();
        final double precursorMZ = inp.getPrecursorMz();
        final double precursorIntensity = 1;
        final String id = inp.getId();
        String peptide = null; // todo do better
        final List<IPeak> peaklist = inp.getPeaks();
        List<Peak> newPeaks = peaklist.stream().map(Adapters::fromPeak).collect(Collectors.toList());

        ClusteringSpectrum ret = new ClusteringSpectrum(
                id,
                precursorMZ,
                precursorIntensity,
                (int) precursorCharge,
                newPeaks);

              ret.setPeptide(inp.getProperty(KnownProperties.IDENTIFIED_PEPTIDE_KEY));
         return ret;
    }


    /**
     * convert old peak to new  peak
     *
     * @param pk !null peak
     * @return !null equivalent peak
     */
    public static Peak fromPeak(IPeak pk) {
        return new Peak(pk.getMz(), pk.getIntensity(), pk.getCount());
    }


    /**
     * operations that merge peaks in lists had better preserve the total count
     * useful in testing and debuggine
     *
     * @param lst !null list of peaks
     * @return total count from all peaks
     */
    public static int getTotalCount(List<Peak> lst) {
        int total = lst.stream().filter(Objects::nonNull).mapToInt(Peak::getCount).sum();
        return total;
    }


}

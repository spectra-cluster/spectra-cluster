package uk.ac.ebi.pride.spectracluster.consensus;

import uk.ac.ebi.pride.spectracluster.io.ParserUtilities;
import uk.ac.ebi.pride.spectracluster.similarity.AllPeaksDotProduct;
import uk.ac.ebi.pride.spectracluster.spectrum.IPeak;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.spectrum.Spectrum;
import uk.ac.ebi.pride.spectracluster.util.Defaults;
import uk.ac.ebi.pride.spectracluster.util.function.peak.BinnedHighestNPeakFunction;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

/**
 * @author Rui Wang
 * @version $Id$
 */
public class ConsensusSpectrumBuilderComparisonMain {

    public static final double DEFAULT_COMPARISON_THRESHOLD = 0.95;

    public static void main(String[] args) {
        if (args.length < 1) {
            System.err.println("Usage: [Input folder]");
            System.exit(1);
        }

        // get all the files from input folder
        File folder = new File(args[0]);
        if (!folder.isDirectory() && !folder.exists()) {
            System.err.println("Input folder must be an existing directory: " + folder.getAbsolutePath());
            System.exit(1);
        }

        // iterate over each file
        for (File file : folder.listFiles()) {
            if (file.getName().endsWith("mgf")) {
                compareConsensusSpectrumBuilder(file);
            }
        }

    }

    private static void compareConsensusSpectrumBuilder(File spectrumFile) {
        // filter the spectra first
        BinnedHighestNPeakFunction filter = new BinnedHighestNPeakFunction();

        // load and filter the spectra
        ISpectrum[] spectra= ParserUtilities.readMGFScans(spectrumFile);
        List<ISpectrum> filteredSpectra = new ArrayList<ISpectrum>(spectra.length);

        for (ISpectrum s : spectra) {
            List<IPeak> peaks = s.getPeaks();
            List<IPeak> filteredPeaks = filter.apply(peaks);

            Spectrum filteredSpectrum = new Spectrum(s, filteredPeaks);

            filteredSpectra.add(filteredSpectrum);
        }

        // create the consensus spectra
        final IConsensusSpectrumBuilder builder = Defaults.getDefaultConsensusSpectrumBuilder();

        builder.addSpectra(spectra);
        ISpectrum unfilteredConsensusSpectrum = builder.getConsensusSpectrum();

        builder.clear();
        builder.addSpectra(filteredSpectra.toArray(new ISpectrum[filteredSpectra.size()]));
        ISpectrum filteredConsensusSpectrum = builder.getConsensusSpectrum();

        // compare the spectra
        final AllPeaksDotProduct allPeaksDotProduct = new AllPeaksDotProduct(Defaults.getSimilarityMZRange());

        final double similarity = allPeaksDotProduct.assessSimilarity(unfilteredConsensusSpectrum, filteredConsensusSpectrum);


        if (similarity < DEFAULT_COMPARISON_THRESHOLD) {
            System.out.println("Consensus spectrum builder comparison failed for file : " + spectrumFile.getAbsolutePath() + " with similarity: " + similarity);
        }
    }
}

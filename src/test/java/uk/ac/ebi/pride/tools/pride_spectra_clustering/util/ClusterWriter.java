package uk.ac.ebi.pride.tools.pride_spectra_clustering.util;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

/**
 * uk.ac.ebi.pride.tools.pride_spectra_clustering.util.ClusterWriter
 *
 * @author Steve Lewis
 */
public class ClusterWriter {
    public static ClusterWriter[] EMPTY_ARRAY = {};
    public static Class THIS_CLASS = ClusterWriter.class;

    /**
     * comparator to compare peaks by intensity first then mz rather than the
     * standard mz then intensity
     */
    public static final Comparator<Peak> BY_MZ = (o1, o2) -> {
        if (o2 == o1)
            return 0;
        if (o1.getMz() != o2.getMz()) {
            return o1.getMz() < o2.getMz() ? -1 : 1;
        }
        if (o1.getIntensity() != o2.getIntensity()) {
            return o2.getIntensity() < o1.getIntensity() ? -1 : 1;
        }
        if (o1.getCount() != o2.getCount()) {
            return o2.getCount() < o1.getCount() ? -1 : 1;
        }

        return 0;
    };


    /**
     * comparator to Spectra by charge then mz
     * standard mz then id
     */
    public static final Comparator<ClusteringSpectrum> BY_SPECTRUM_MZ = (o1, o2) -> {
        if (o2 == o1)
            return 0;
        int charge1 = (int) (0.5 + o1.getPrecursorCharge());
        int charge2 = (int) (0.5 + o2.getPrecursorCharge());
        if (charge1 != charge2) {
            return charge1 < charge2 ? -1 : 1;
        }
        if (o1.getPrecursorMZ() != o2.getPrecursorMZ()) {
            return o1.getPrecursorMZ() < o2.getPrecursorMZ() ? -1 : 1;
        }
        return o1.getId().compareTo(o2.getId());
    };

    /**
     * comparator to SpectraCluster by charge then mz
     * then intensity
     */
    public static final Comparator<SpectraCluster> BY_CHARGE_THEN_MZ = (o1, o2) -> {
        if (o2 == o1)
            return 0;
        int charge1 = (int) (0.5 + o1.getAverageCharge());
        int charge2 = (int) (0.5 + o2.getAverageCharge());
        if (charge1 != charge2) {
            return charge1 < charge2 ? -1 : 1;
        }
        if (o1.getAverageMz() != o2.getAverageMz()) {
            return o1.getAverageMz() < o2.getAverageMz() ? -1 : 1;
        }
        if (o1.getAverageIntensity() != o2.getAverageIntensity()) {
            return o2.getAverageIntensity() < o1.getAverageIntensity() ? -1 : 1;
        }

        return 0;
    };


    public static void dumpSpectra(String outFile, List<SpectraCluster> generatedClusterX) {
        try {
            // defensive copy
            List<SpectraCluster> generatedCluster = new ArrayList<>(generatedClusterX);
            generatedCluster.sort(BY_CHARGE_THEN_MZ);

            PrintWriter out = new PrintWriter(new FileWriter(outFile));
            for (SpectraCluster spectraCluster : generatedCluster) {
                dumpCluster(spectraCluster, out);
            }
            out.close();
        } catch (IOException e) {
            throw new RuntimeException(e);
        }

    }

    private static int gClusterId = 1000;

    public static void dumpCluster(SpectraCluster spectraCluster, PrintWriter out) {
        List<ClusteringSpectrum> spectrumsX = spectraCluster.getSpectra();
        if (spectrumsX.size() == 5)
            spectrumsX = spectraCluster.getSpectra();
        String id = Integer.toString(gClusterId++);
        int charge = (int) (0.5 + spectraCluster.getAverageCharge());
        double mz = spectraCluster.getAverageMz();
        out.append("BEGIN CLUSTER");
        out.append(" Id=").append(id);
        out.append(" Charge=").append(String.valueOf(charge));
        out.append("\n");
        List<Peak> consensusSpectrum = spectraCluster.consensusSpectrum;

        dumpSpectrum(id, mz, charge, consensusSpectrum, out);

        // defensive copy
        List<ClusteringSpectrum> spectrums = new ArrayList<>(spectrumsX);
        spectrums.sort(BY_SPECTRUM_MZ);

        for (ClusteringSpectrum spectrum : spectrums) {
            dumpSpectrum(spectrum, out);
        }

        out.append("END CLUSTER");
        out.append("\n");


    }


    public static void dumpSpectrum(String id, double mz, double charge, List<Peak> peaklistX, PrintWriter out) {
        out.append("BEGIN IONS");
        out.append("\n");

        out.append("TITLE=").append(id);
        out.append("\n");


        out.append("PEPMASS=").append(String.valueOf(mz));
        out.append("\n");

        out.append("CHARGE=").append(String.valueOf(charge));
        if (charge > 0)
            out.append("+");
        out.append("\n");

        dumpPeaks(peaklistX, out);


        out.append("END IONS");
        out.append("\n");


    }


    public static void dumpSpectrum(ClusteringSpectrum sp, PrintWriter out) {
        out.append("BEGIN IONS");
        out.append("\n");

        out.append("TITLE=").append(sp.getId());
        out.append("\n");


        out.append("PEPMASS=").append(String.valueOf(sp.getPrecursorMZ()));
        out.append("\n");

        Integer precursorCharge = sp.getPrecursorCharge();
        out.append("CHARGE=").append(String.valueOf(precursorCharge));
        if (precursorCharge > 0)
            out.append("+");
        out.append("\n");

        List<Peak> peaklistX = sp.getPeaklist();
        dumpPeaks(peaklistX, out);


        out.append("END IONS");
        out.append("\n");


    }

    public static void dumpPeaks(List<Peak> peaklistX, PrintWriter out) {
        // defensive copy
        List<Peak> peaklist = new ArrayList<>(peaklistX);
        peaklist.sort(BY_MZ);

        for (Peak peak : peaklist) {
            String item = String.format("%10.5f", peak.getMz()).trim();
            out.print(item);
            out.print("\t");
            String item2 = String.format("%8.2f", peak.getIntensity()).trim();
            out.print(item2);
            out.print("\t");
            String item3 = Integer.toString(peak.getCount());
            out.print(item3);
            out.print("\n");

        }
    }


    public static void dumpPeaks(List<Peak> peaklistX, PrintStream out) {
        // defensive copy
        List<Peak> peaklist = new ArrayList<>(peaklistX);
        peaklist.sort(BY_MZ);

        for (Peak peak : peaklist) {
            String item = String.format("%10.5f", peak.getMz()).trim();
            out.print(item);
            out.print("\t");
            String item2 = String.format("%8.2f", peak.getIntensity()).trim();
            out.print(item2);
            out.print("\t");
            String item3 = Integer.toString(peak.getCount());
            out.print(item3);
            out.print("\n");

        }
    }

}

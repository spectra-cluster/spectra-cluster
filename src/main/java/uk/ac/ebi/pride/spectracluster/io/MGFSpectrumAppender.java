package uk.ac.ebi.pride.spectracluster.io;

import uk.ac.ebi.pride.spectracluster.spectrum.IPeak;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.spectrum.KnownProperties;

import java.io.IOException;
import java.util.Properties;

/**
 * Write spectrum out as MGF peak list format
 *
 * @author Rui Wang
 * @version $Id$
 */
public class MGFSpectrumAppender implements ISpectrumAppender {
    public static MGFSpectrumAppender INSTANCE = new MGFSpectrumAppender();

    protected MGFSpectrumAppender() {
    }

    @Override
    public void appendSpectrum(Appendable out, ISpectrum spectrum, Object... otherData) {
        try {
            out.append("BEGIN IONS");
            out.append("\n");

            appendTitle(spectrum, out);
            out.append("\n");

            double precursorCharge = spectrum.getPrecursorCharge();
            double massChargeRatio = spectrum.getPrecursorMz();

            out.append("PEPMASS=").append(String.valueOf(massChargeRatio));
            out.append("\n");

            out.append("CHARGE=").append(String.valueOf(precursorCharge));
            if (precursorCharge > 0)
                out.append("+");
            out.append("\n");

            appendProperties(spectrum, out);

            appendPeaks(spectrum, out);

            out.append("END IONS");
            out.append("\n");
        } catch (IOException e) {
            throw new AppenderException(e);
        }
    }

    public void appendProperties(ISpectrum spectrum, Appendable out) {
        final Properties properties = spectrum.getProperties();
        try {
            for (String s : properties.stringPropertyNames()) {
                final String property = properties.getProperty(s);
                final String line = KnownProperties.toMGFLine(s, property);
                out.append(line);
                out.append("\n");
            }
        } catch (IOException e) {
            throw new UnsupportedOperationException(e);
        }
    }

    /**
     * override to add peptide later
     *
     * @param out
     * @throws IOException
     */
    public void appendTitle(final ISpectrum spectrum, final Appendable out) throws IOException {
        out.append("TITLE=id=").append(spectrum.getId());
//        final String peptide = spectrum.getProperty(KnownProperties.IDENTIFIED_PEPTIDE_KEY);
//        if (peptide != null && peptide.length() > 0)
//            out.append(",sequence=").append(peptide);
    }

    protected void appendPeaks(final ISpectrum spectrum, final Appendable out) throws IOException {
        for (IPeak peak : spectrum.getPeaks()) {
            String line = String.format("%10.3f", peak.getMz()).trim() + "\t" +
                    String.format("%10.3f", peak.getIntensity()).trim();
            out.append(line);
            out.append("\n");
        }
    }

    @Override
    public void appendStart(Appendable out, Object... otherData) {

    }

    @Override
    public void appendEnd(Appendable out, Object... otherData) {

    }
}

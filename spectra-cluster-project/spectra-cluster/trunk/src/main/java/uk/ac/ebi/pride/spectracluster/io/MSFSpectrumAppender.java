package uk.ac.ebi.pride.spectracluster.io;

import uk.ac.ebi.pride.spectracluster.spectrum.IPeak;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;

import java.io.IOException;

/**
 * Spectrum appender for MSF file format
 *
 * @author Rui Wang
 * @version $Id$
 */
public class MSFSpectrumAppender implements ISpectrumAppender {

    public static MSFSpectrumAppender INSTANCE = new MSFSpectrumAppender();

    private MSFSpectrumAppender() {}

    @Override
    public void appendSpectrum(Appendable out, ISpectrum spectrum, Object... otherData) {
        try {
            appendMSFComment(spectrum, out);
            appendMSFPeaks(spectrum, out);
            out.append("\n");
        } catch (IOException e) {
            throw new AppenderException(e);

        }
    }

    private void appendMSFComment(final ISpectrum spectrum, final Appendable out) {
        //todo: the content of this method is Steve's original implementation, need to find an alternative
//        List<String> strings = new ArrayList<String>(properties.keySet());
//        Collections.sort(strings);
//        boolean first = true;
//        try {
//            out.append(ParserUtilities.COMMENT_START + " ");
//            for (String string : strings) {
//                if (!first) {
//                    out.append(" ");
//                    first = false;
//                } else {
//                    first = false;
//                }
//                out.append(string);
//                out.append("=");
//                String value = properties.get(string);
//                if (value.contains(" ")) {
//                    value = "\"" + value + "\"";
//                }
//                out.append(value);
//            }
//            out.append("\n");
//        } catch (IOException e) {
//            throw new RuntimeException(e);
//
//        }
    }

    private void appendMSFPeaks(final ISpectrum spectrum, final Appendable out) throws IOException {
        out.append("Num peaks: ").append(String.valueOf(spectrum.getPeaks().size())).append("\n");
        appendPeaks(spectrum, out);
    }

    private void appendPeaks(final ISpectrum spectrum, final Appendable out) throws IOException {
        for (IPeak peak : spectrum.getPeaks()) {
            String line = String.format("%10.3f", peak.getMz()).trim() + "\t" +
                    String.format("%10.3f", peak.getIntensity()).trim();
            out.append(line);
            out.append("\n");
        }
    }

    @Override
    public void appendStart(Appendable out, Object... otherData) {
        //todo: the content of this method is Steve's original implementation, need to find an alternative
//        String peptide = getPeptide();
//        try {
//            out.append(ParserUtilities.NAME_START + " " + peptide + "/" + getPrecursorCharge());
//            out.append("\n");
//            String id = getId();
//            if (id != null) {
//                out.append(ParserUtilities.LIBID_START + " " + getId());
//                out.append("\n");
//
//            }
//            //        out.append(ParserUtilities.MW_START + String.format("%10.3f", getPrecursorMz() * getPrecursorCharge()).trim());
//            //        out.append("\n");
//            out.append(ParserUtilities.PRECURSORMZ_START + " " + String.format("%10.3f", getPrecursorMz()).trim());
//            out.append("\n");
//
//            String prop = getProperty("Status");
//            if (prop != null) {
//                out.append(ParserUtilities.STATUS_START + prop);
//                out.append("\n");
//            }
//            prop = getProperty("molecularWeight");
//            if (prop != null) {
//                out.append(ParserUtilities.MW_START + prop);
//                out.append("\n");
//            }
//            prop = getProperty("FullName");
//            if (prop != null) {
//                out.append(ParserUtilities.FULL_NAME_START + prop);
//                out.append("\n");
//            }
//        } catch (IOException e) {
//            throw new RuntimeException(e);
//
//        }
    }

    @Override
    public void appendEnd(Appendable out, Object... otherData) {

    }
}

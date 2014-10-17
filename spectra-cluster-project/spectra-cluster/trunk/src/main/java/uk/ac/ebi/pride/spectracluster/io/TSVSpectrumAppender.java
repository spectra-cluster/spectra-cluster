package uk.ac.ebi.pride.spectracluster.io;

import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;

import java.io.IOException;

/**
 * Append spectrum as TSV file
 *
 * @author Rui Wang
 * @version $Id$
 */
public class TSVSpectrumAppender implements ISpectrumAppender {

    @Override
    public void appendSpectrum(Appendable out, ISpectrum spectrum, Object... otherData) throws AppenderException {
        try {
            out.append(spectrum.getId());
            out.append("\t");
            out.append(Integer.toString(spectrum.getPrecursorCharge()));
            out.append("\t");
            String mzString = String.format("%10.2f", spectrum.getPrecursorMz()).trim();
            out.append(mzString);
            out.append("\n");
        } catch (IOException e) {
            throw new AppenderException(e);

        }
    }

    @Override
    public void appendStart(Appendable out, Object... otherData) throws AppenderException {

    }

    @Override
    public void appendEnd(Appendable out, Object... otherData) throws AppenderException {

    }
}

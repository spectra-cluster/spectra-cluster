package uk.ac.ebi.pride.spectracluster.io;

import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;

import java.io.IOException;

/**
 * This code is licensed under the Apache License, Version 2.0 (the
 * "License"); you may not use this file except in compliance
 * with the License.  You may obtain a copy of the License at
 * <p>
 * http://www.apache.org/licenses/LICENSE-2.0
 * <p>
 * ==Overview==
 * <p>
 * This class
 * <p>
 * Created by ypriverol (ypriverol@gmail.com) on 08/12/2017.
 */
public class MGFSpectrumAppenderIntCharge extends MGFSpectrumAppender{

    public static MGFSpectrumAppenderIntCharge INSTANCE = new MGFSpectrumAppenderIntCharge();

    protected MGFSpectrumAppenderIntCharge() {
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

            out.append("CHARGE=").append(String.valueOf((int) precursorCharge));
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
}

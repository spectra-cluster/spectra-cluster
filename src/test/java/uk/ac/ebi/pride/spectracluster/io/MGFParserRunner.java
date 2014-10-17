package uk.ac.ebi.pride.spectracluster.io;

import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.LineNumberReader;
import java.util.ArrayList;
import java.util.List;

/**
 * @author Rui Wang
 * @version $Id$
 */
public class MGFParserRunner {
    public static void main(String[] args) throws IOException {

        LineNumberReader reader = null;
        try {
            reader = new LineNumberReader(new FileReader(new File(args[0])));
            List<ISpectrum> spectra = new ArrayList<ISpectrum>();
            ISpectrum spectrum = null;
            while((spectrum = ParserUtilities.readMGFScan(reader)) != null) {
                System.out.println(spectrum.getId());
                spectra.add(spectrum);
            }


        } finally {
            if (reader != null)
                reader.close();
        }
    }
}

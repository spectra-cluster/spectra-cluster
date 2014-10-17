package uk.ac.ebi.pride.spectracluster.util;


import uk.ac.ebi.pride.tools.jmzreader.*;
import uk.ac.ebi.pride.tools.jmzreader.model.*;
import uk.ac.ebi.pride.tools.jmzreader.model.Spectrum;
import uk.ac.ebi.pride.tools.mgf_parser.*;
import uk.ac.ebi.pride.tools.mgf_parser.model.*;

import java.io.*;
import java.net.*;
import java.util.*;

/**
 * uk.ac.ebi.pride.spectracluster.util.JMZTabUtilities
 *
 * @author Steve Lewis
 * @date 29/05/2014
 */
public class JMZTabUtilities {

    /**
     * read a resource mgf as a list of spectra
     *
     * @param resName
     * @return
     */
    public static List<Spectrum> readSpectrumsFromResource(String resName) {
        try {
            File specFile;
            MgfFile mgfFile;
            List<Spectrum> spectra;
            URL testFile = ClusteringTestUtilities.class.getClassLoader().getResource(resName);

            assert testFile != null;
            specFile = new File(testFile.toURI());

            mgfFile = new MgfFile(specFile);

            spectra = new ArrayList<Spectrum>(mgfFile.getMs2QueryCount());
            Iterator<Ms2Query> it = mgfFile.getMs2QueryIterator();
            while (it.hasNext()) {
                Ms2Query query = it.next();
                if (query.getPrecursorIntensity() == null)
                    query.setPeptideIntensity(1.0);

                spectra.add(query);
            }
            return spectra;
        } catch (URISyntaxException e) {
            throw new RuntimeException(e);
        } catch (JMzReaderException e) {
            throw new RuntimeException(e);
        }


    }

    /**
     * read a resource mgf as a list of spectra
     *
     * @return
     */
    public static List<Spectrum> readSpectrumsFromResource() {
        return readSpectrumsFromResource(ClusteringTestUtilities.SAMPLE_MGF_FILE);
    }
}

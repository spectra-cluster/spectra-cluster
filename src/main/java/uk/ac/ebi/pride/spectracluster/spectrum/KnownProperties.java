package uk.ac.ebi.pride.spectracluster.spectrum;

import java.util.Collections;
import java.util.HashMap;
import java.util.Map;
import java.util.Properties;

/**
 * uk.ac.ebi.pride.spectracluster.spectrum.KnownProperties
 * This class has properties of spectra and clusters known well enough to be hard coded and hints about
 * how they are written to MGF and CGF files
 *
 * see http://www.matrixscience.com/help/data_file_help.html
 * for defined MGF Keys
 *
 * @author Johannes Griss
 * @author Steve Lewis
 * @author Yasset Perez-Riverol
 */
public class KnownProperties {

    // Known properties keys
    public static final String IDENTIFIED_PEPTIDE_KEY = "identifiedPeptide";
    public static final String ANNOTATION_KEY = "annotation";
    public static final String TAXONOMY_KEY = "accession";
    public static final String PROTEIN_KEY = "protein"; // database: protein
    public static final String MODIFICATION_KEY = "modification"; // database: protein
    public static final String INSTRUMENT_KEY = "instrument";
    public static final String SPECTRUM_TITLE = "custom_title";
    public static final String RETENTION_TIME = "retention_time";
    public static final String PSM_DECOY_STATUS = "psm_decoy_status";
    public static final String PSM_FDR_SCORES = "psm_fdr_scores";
    public static final String MIN_COMPARISONS = "min_comp";
    public final static String ADDING_SCORE = "adding_score";

    public static final String IDENTIFIED_PEPTIDE_MGF_KEY = "SEQ";
    public static final String ANNOTATION_MGF_KEY = "USER00";
    public static final String TAXONOMY_MGF_KEY = "TAXONOMY";
    public static final String PROTEIN_MGF_KEY = "USER02";
    public static final String MODIFICATION_MGF_KEY = "USER03";
    public static final String INSTRUMENT_MGF_KEY = "INSTRUMENT";
    public static final String SPECTRUM_MGF_TITLE = "USER04";
    public static final String DECOY_STATUS_MGF_KEY = "USER05";
    public static final String FDR_SCORES_MGF_KEY = "USER06";
    public static final String RETENTION_TIME_MGF_KEY = "RTINSECONDS";
    public static final String MIN_COMPARISONS_MGF_KEY = "MIN_COMP";
    public static final String ADDING_SCORE_MGF_KEY = "ADDING_SCORE";


    public static final String UNKNOWN_MGF_KEY = "USER12";

    // ===========================
    // Known cluster Properties
    public static final String MOST_COMMON_PEPTIDE_KEY = "mostCommonPeptide";
    public static final String N_HIGHEST_PEAKS = "n_highest_peaks"; // highest peaks to use for comparison


    // =====================
    /**
     * this section related to tags in MGF files where
     * SEQ, USER00, USER01, USER02 .. User12 are allowed
     */
    private static Map<String, String> INTERNAL_KEY_TO_MGF_KEY = new HashMap<String, String>();
    private static Map<String, String> INTERNAL_MGF_KEY_TO_KEY = new HashMap<String, String>();

    static {
        INTERNAL_KEY_TO_MGF_KEY.put(IDENTIFIED_PEPTIDE_KEY, IDENTIFIED_PEPTIDE_MGF_KEY);
        INTERNAL_KEY_TO_MGF_KEY.put(ANNOTATION_KEY, ANNOTATION_MGF_KEY);
        INTERNAL_KEY_TO_MGF_KEY.put(TAXONOMY_KEY, TAXONOMY_MGF_KEY);
        INTERNAL_KEY_TO_MGF_KEY.put(PROTEIN_KEY, PROTEIN_MGF_KEY);
        INTERNAL_KEY_TO_MGF_KEY.put(MODIFICATION_KEY, MODIFICATION_MGF_KEY);
        INTERNAL_KEY_TO_MGF_KEY.put(INSTRUMENT_KEY, INSTRUMENT_MGF_KEY);
        INTERNAL_KEY_TO_MGF_KEY.put(SPECTRUM_TITLE, SPECTRUM_MGF_TITLE);
        INTERNAL_KEY_TO_MGF_KEY.put(RETENTION_TIME, RETENTION_TIME_MGF_KEY);
        INTERNAL_KEY_TO_MGF_KEY.put(PSM_DECOY_STATUS, DECOY_STATUS_MGF_KEY);
        INTERNAL_KEY_TO_MGF_KEY.put(PSM_FDR_SCORES, FDR_SCORES_MGF_KEY);
        INTERNAL_KEY_TO_MGF_KEY.put(MIN_COMPARISONS, MIN_COMPARISONS_MGF_KEY);
        INTERNAL_KEY_TO_MGF_KEY.put(ADDING_SCORE, ADDING_SCORE_MGF_KEY);

        INTERNAL_MGF_KEY_TO_KEY.put(IDENTIFIED_PEPTIDE_MGF_KEY, IDENTIFIED_PEPTIDE_KEY);
        INTERNAL_MGF_KEY_TO_KEY.put(ANNOTATION_MGF_KEY, ANNOTATION_KEY);
        INTERNAL_MGF_KEY_TO_KEY.put(TAXONOMY_MGF_KEY, TAXONOMY_KEY);
        INTERNAL_MGF_KEY_TO_KEY.put(PROTEIN_MGF_KEY, PROTEIN_KEY);
        INTERNAL_MGF_KEY_TO_KEY.put(MODIFICATION_MGF_KEY, MODIFICATION_KEY);
        INTERNAL_MGF_KEY_TO_KEY.put(INSTRUMENT_MGF_KEY, INSTRUMENT_KEY);
        INTERNAL_MGF_KEY_TO_KEY.put(SPECTRUM_MGF_TITLE, SPECTRUM_TITLE);
        INTERNAL_MGF_KEY_TO_KEY.put(RETENTION_TIME_MGF_KEY, RETENTION_TIME);
        INTERNAL_MGF_KEY_TO_KEY.put(DECOY_STATUS_MGF_KEY, PSM_DECOY_STATUS);
        INTERNAL_MGF_KEY_TO_KEY.put(FDR_SCORES_MGF_KEY, PSM_FDR_SCORES);
        INTERNAL_MGF_KEY_TO_KEY.put(MIN_COMPARISONS_MGF_KEY, MIN_COMPARISONS);
        INTERNAL_MGF_KEY_TO_KEY.put(ADDING_SCORE_MGF_KEY, ADDING_SCORE);

    }

    public static Map<String, String> KEY_TO_MGF_KEY = Collections.unmodifiableMap(INTERNAL_KEY_TO_MGF_KEY);
    public static Map<String, String> MGF_KEY_TO_KEY = Collections.unmodifiableMap(INTERNAL_MGF_KEY_TO_KEY);

    /**
     * take property - value pair to a line to insert in MGF
     * USER12 is any unknown
     *
     * @param property Property name
     * @param value The property's value
     * @return The encoded property / value pair based on the names of known properties.
     */
    public static String toMGFLine(String property, String value) {
        String key = KEY_TO_MGF_KEY.get(property);
        if (key == null) {
            return UNKNOWN_MGF_KEY + "=" + property + "=" + value;
        } else {
            return key + "=" + value;
        }
    }

    /**
     * parse an mgf line like SEQ or USER00..USER12
     * @param props  properties to add
     * @param line The MGF line to test.
     * @return   true if successfully handled
     */
    public static boolean addMGFProperties(Properties props,String line) {
        if(line.contains("="))  {
            int index = line.indexOf("=") ;
            String key = line.substring(0,index);
            String value = line.substring(index + 1,line.length());
            return handleKnownProperty(props,key, value);
        }
        else {
            return false;
        }
    }

    private static boolean handleUnknownProperty(Properties props,String key1, String key2, String value) {
          if (!UNKNOWN_MGF_KEY.equals(key1))
              return false;
       //     throw new UnsupportedOperationException("Properties need to be  USER12=name=value not " + key1);
        props.setProperty(key2,value);
        return true;
    }

    private static boolean handleKnownProperty(Properties props,String key, String value) {
        if (key.startsWith("USER12="))
            key = key.substring(7);

        String realKey =  MGF_KEY_TO_KEY.get(key);
        if (realKey == null)
            return false;
       //      throw new UnsupportedOperationException("Properties need to be known key= value" + key);
        props.setProperty(realKey,value);
        return true;
    }
}

package uk.ac.ebi.pride.spectracluster.spectrum;

import uk.ac.ebi.pride.spectracluster.util.*;

import java.io.*;

/**
 * IPeak is an interface which represents a peak in a spectrum
 *
 * @author Johannes Griss
 * @author Steve Lewis
 * @author Rui Wang
 */
public interface IPeak extends Equivalent<IPeak>, Comparable<IPeak>, Serializable {

    public static final IPeak[] EMPTY_ARRAY = {};

    /**
     * Peak m/z
     */
    public float getMz();

    /**
     * Peak intensity
     */
    public float getIntensity();

    /**
     * If the peak is part of a consensus spectrum this number represents the number of
     * spectra making up the consensus spectrum that contain the respective peak. In normal spectra
     * this number is always 1.
     */
    public int getCount();


}

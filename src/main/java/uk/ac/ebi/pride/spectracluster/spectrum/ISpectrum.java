package uk.ac.ebi.pride.spectracluster.spectrum;

import uk.ac.ebi.pride.spectracluster.util.Equivalent;

import java.util.List;
import java.util.Properties;

/**
 * uk.ac.ebi.pride.spectracluster.spectrum.IPeaksSpectrum
 * used by Spectra like get highest peaks which are incomplete
 * spectra
 *
 * @author Steve Lewis
 * @author Rui Wang
 *         Date: 6/20/13
 */
public interface ISpectrum extends ISpectrumQuality, Equivalent<ISpectrum>, Comparable<ISpectrum> {

    /**
     * globally unique id
     *
     * @return !null id
     */
    String getId();

    /**
     * get precursor m/z
     */
    float getPrecursorMz();

    /**
     * get charge - mixed charge
     */
    int getPrecursorCharge();

    /**
     * return the sum of all intensities
     */
    double getTotalIntensity();

    /**
     * return the sum  Square of all intensities
     */
    double getSumSquareIntensity();

    /**
     * return unmodifiable peaks sorted by MZ
     *
     * @return !null array of peaks
     */
    List<IPeak> getPeaks();

    /**
     * return number of peaks
     *
     * @return count
     */
    int getPeaksCount();

    /**
     * get the highest intensity peaks sorted by MZ - this value may be cached
     *
     * @param numberRequested number peaks requested
     * @return Peaks spectrum
     */
    ISpectrum getHighestNPeaks(int numberRequested);

    /**
     * return as a spectrum the highest n peaks as defined in majorPeakCount
     * this follows Frank et all suggestion that all spectra in a cluster will share at least one of these
     *
     * @param majorPeakCount The number of highest peaks to consider "major"
     * @return
     */
    int[] asMajorPeakMZs(int majorPeakCount);

    /**
     * does the spectrum contain this is a major peak
     *
     * @param mz             peak as int
     * @param majorPeakCount The number of highest peaks to consider "major"
     * @return true if so
     */
    boolean containsMajorPeak(final int mz, int majorPeakCount);

    /**
     * return a property of null if none exists
     * look in ISpectrum for known keys
     *
     * @param key
     * @return
     */
    String getProperty(String key);

    /**
     * look in ISpectrum for known keys
     *
     * @param key
     * @param value
     */
    void setProperty(String key, String value);

    /**
     * Only for internal use in copy constructor
     * Note this is not safe
     * This is not really deprecated but it warns only for
     * internal use
     */
    @Deprecated
    Properties getProperties();


}
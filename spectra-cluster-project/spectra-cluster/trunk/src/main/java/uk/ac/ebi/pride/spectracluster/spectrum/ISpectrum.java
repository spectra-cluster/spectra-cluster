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
    public String getId();

    /**
     * get precursor m/z
     */
    public float getPrecursorMz();

    /**
     * get charge - mixed charge
     */
    public int getPrecursorCharge();

    /**
     * return the sum of all intensities
     */
    public double getTotalIntensity();

    /**
     * return the sum  Square of all intensities
     */
    public double getSumSquareIntensity();

    /**
     * return unmodifiable peaks sorted by MZ
     *
     * @return !null array of peaks
     */
    public List<IPeak> getPeaks();

    /**
     * return number of peaks
     *
     * @return count
     */
    public int getPeaksCount();

    /**
     * get the highest intensity peaks sorted by MZ - this value may be cached
     *
     * @param numberRequested number peaks requested
     * @return Peaks spectrum
     */
    public ISpectrum getHighestNPeaks(int numberRequested);

    /**
     * return as a spectrum the highest n peaks as defined in majorPeakCount
     * this follows Frank et all suggestion that all spectra in a cluster will share at least one of these
     *
     * @param majorPeakCount The number of highest peaks to consider "major"
     * @return
     */
    public int[] asMajorPeakMZs(int majorPeakCount);

    /**
     * does the spectrum contain this is a major peak
     *
     * @param mz             peak as int
     * @param majorPeakCount The number of highest peaks to consider "major"
     * @return true if so
     */
    public boolean containsMajorPeak(final int mz, int majorPeakCount);

    /**
     * return a property of null if none exists
     * look in ISpectrum for known keys
     *
     * @param key
     * @return
     */
    public String getProperty(String key);

    /**
     * look in ISpectrum for known keys
     *
     * @param key
     * @param value
     */
    public void setProperty(String key, String value);

    /**
     * Only for internal use in copy constructor
     * Note this is not safe
     * This is not really deprecated but it warns only for
     * internal use
     */
    @Deprecated
    public Properties getProperties();


}
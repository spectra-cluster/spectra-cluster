package uk.ac.ebi.pride.spectracluster.cluster;

import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.util.Equivalent;

import javax.annotation.Nonnull;
import java.io.*;
import java.util.*;

/**
 * look in ClusterUtilities and ClusterSimilarityUtilities for older methods
 * @author Rui Wang
 * @version $Id$
 */
public interface ICluster extends ISpectrumHolder,
        Equivalent<ICluster>,
        Comparable<ICluster>, Serializable {

    /**
     * Get cluster id
     */
    public String getId();

    /**
     * build an id from spectral ids
     *
     * @return
     */
    public String getSpectralId();

    /**
     * concensus spectrum MZ. If not available (ie. no spectra in cluster)
     * 0 is returned.
     *
     * @return
     */
    public float getPrecursorMz();

    /**
     * concensus spectrum Charge. If not available (ie. no spectra in cluster)
     * 0 is returned.
     *
     * @return
     */
    public int getPrecursorCharge();

    /**
     * Get consensus spectrum
     */
    public ISpectrum getConsensusSpectrum();

    /**
     * real spectrum with the highest quality - this is a
     * good way to compare clusters
     *
     * @return !null spectrum
     */
    public ISpectrum getHighestQualitySpectrum();

    /**
     * all internally spectrum
     */
    @Nonnull
    public List<ISpectrum> getHighestQualitySpectra();

    /**
     * all internally spectrum
     */
    @Nonnull
    public List<ISpectrum> getClusteredSpectra();


    /**
     * count of internal spectrum
     */
    public int getClusteredSpectraCount();


    /**
     * return a set of all ids
     *
     * @return
     */
    @Nonnull
    public Set<String> getSpectralIds();

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

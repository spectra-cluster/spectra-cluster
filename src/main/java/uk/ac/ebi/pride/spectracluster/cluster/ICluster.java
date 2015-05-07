package uk.ac.ebi.pride.spectracluster.cluster;

import uk.ac.ebi.pride.spectracluster.consensus.IConsensusSpectrumBuilder;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.util.ComparisonMatch;
import uk.ac.ebi.pride.spectracluster.util.Equivalent;

import javax.annotation.Nonnull;
import java.io.Serializable;
import java.util.List;
import java.util.Properties;
import java.util.Set;

/**
 * look in ClusterUtilities and ClusterSimilarityUtilities for older methods
 *
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
     * Set cluster id
     *
     * @return
     */
    public void setId(String id);

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
     * Get consensus spectrum builder
     * @return  consensus spectrum builder
     */
    public IConsensusSpectrumBuilder getConsensusSpectrumBuilder();

    /**
     * real spectrum with the highest quality - this is a
     * good way to compare clusters
     *
     * @return !null spectrum
     */
    public ISpectrum getHighestQualitySpectrum();

    /**
     * The N highest quality spectra. N depends on the constant set in SpectralQualityHolder
     */
    @Nonnull
    @Deprecated // TODO jg: getHighestQualitySpectra this function does not seem to be used
    public List<ISpectrum> getHighestQualitySpectra(); // TODO jg: getHighestQualitySpectra - check how this function is used

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

    /**
     * Indicates whether the cluster implementation stores peak lists.
     * @return
     */
    public boolean storesPeakLists();

    /**
     * The results of the last N comparisons.
     * @return
     */
    public List<ComparisonMatch> getComparisonMatches();
}

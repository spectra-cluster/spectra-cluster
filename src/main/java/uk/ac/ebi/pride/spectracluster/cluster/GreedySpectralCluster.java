package uk.ac.ebi.pride.spectracluster.cluster;

import uk.ac.ebi.pride.spectracluster.consensus.BinnedGreedyConsensusSpectrum;
import uk.ac.ebi.pride.spectracluster.consensus.IConsensusSpectrumBuilder;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.spectrum.Spectrum;
import uk.ac.ebi.pride.spectracluster.util.CompareTo;
import uk.ac.ebi.pride.spectracluster.util.ComparisonMatch;
import uk.ac.ebi.pride.spectracluster.util.MZIntensityUtilities;
import uk.ac.ebi.pride.spectracluster.util.SpectrumUtilities;

import java.util.*;
import java.util.concurrent.CopyOnWriteArrayList;


/**
 * This implementation if ICluster only supports the addition of
 * spectra. It does not keep the actual spectra after they were added but
 * only their ids.
 *
 * @author Johannes Griss
 */
public class GreedySpectralCluster implements ICluster {
    /**
     * The N (defined here) highest comparison matches will be
     * saved
     */
    public static final int SAVED_COMPARISON_MATCHES = 30;
    private final List<ComparisonMatch> bestComparisonMatches = new ArrayList<ComparisonMatch>(SAVED_COMPARISON_MATCHES);
    private float lowestBestComparisonSimilarity = 0;
    private Set<String> bestComparisonMatchIds = null;

    private String id;
    private final Set<String> spectraIds = new HashSet<String>();
    private final Properties properties = new Properties();
    private final List<SpectrumHolderListener> spectrumHolderListeners = new CopyOnWriteArrayList<SpectrumHolderListener>();

    /**
     * Clustered spectra are only stored WITHOUT their peak lists.
     */
    private List<ISpectrum> clusteredSpectra = new ArrayList<ISpectrum>();

    private final BinnedGreedyConsensusSpectrum consensusSpectrumBuilder;

    public GreedySpectralCluster(String id) {
        this.id = id;
        this.consensusSpectrumBuilder = BinnedGreedyConsensusSpectrum.FACTORY.getGreedyConsensusSpectrumBuilder(id);
        addSpectrumHolderListener(this.consensusSpectrumBuilder);
    }

    public GreedySpectralCluster(ICluster cluster) {
        this.id = cluster.getId();

        // copy the basic parameters
        this.properties.putAll(cluster.getProperties());

        // copy the GreedySpectralCluster specific properties
        if (GreedySpectralCluster.class.isInstance(cluster)) {
            GreedySpectralCluster existingCluster = (GreedySpectralCluster) cluster;

            // copy the comparison matches
            this.bestComparisonMatches.addAll(existingCluster.bestComparisonMatches);
            this.lowestBestComparisonSimilarity = existingCluster.lowestBestComparisonSimilarity;
            this.bestComparisonMatchIds = existingCluster.bestComparisonMatchIds;
            // for greedy clusters the consensus spectrum must be copied since it cannot be derived from the actual spectra
            this.consensusSpectrumBuilder = (BinnedGreedyConsensusSpectrum) existingCluster.getConsensusSpectrumBuilder();
            addSpectrumHolderListener(this.consensusSpectrumBuilder);

            this.clusteredSpectra.addAll(existingCluster.getClusteredSpectra()); // peak lists are already removed
            this.spectraIds.clear();
            for (ISpectrum spectrum : clusteredSpectra) {
                spectraIds.add(spectrum.getId());
            }
        } else {
            // rebuild with a GreedyConsensusSpectrum
            this.consensusSpectrumBuilder = BinnedGreedyConsensusSpectrum.FACTORY.getGreedyConsensusSpectrumBuilder(id);
            addSpectrumHolderListener(this.consensusSpectrumBuilder);

            if (!cluster.storesPeakLists())
                throw new IllegalStateException("Cannot copy cluster without peak lists that is not a GreedyCluster.");

            ISpectrum[] existingSpectra = new ISpectrum[cluster.getClusteredSpectraCount()];
            existingSpectra = cluster.getClusteredSpectra().toArray(existingSpectra);
            addSpectra(existingSpectra);
        }
    }

    public GreedySpectralCluster(String id, List<ISpectrum> clusteredSpectra, BinnedGreedyConsensusSpectrum consensusSpectrumBuilder, List<ComparisonMatch> bestComparisonMatches) {
        this.id = id;
        this.clusteredSpectra = clusteredSpectra;
        this.consensusSpectrumBuilder = consensusSpectrumBuilder;

        addSpectrumHolderListener(this.consensusSpectrumBuilder);
        setComparisonMatches(bestComparisonMatches);

        getSpectralIds();
    }

    /**
     * return a set of all ids
     *
     * @return A set of strings representing all Ids.
     */
    @Override
    public Set<String> getSpectralIds() {
        if (this.spectraIds.isEmpty()) {
            List<ISpectrum> clusteredSpectra1 = getClusteredSpectra();
            for (ISpectrum iSpectrum : clusteredSpectra1) {
                spectraIds.add(iSpectrum.getId());
            }
        }
        return Collections.unmodifiableSet(spectraIds);
    }


    /**
     * add a change listener
     * final to make sure this is not duplicated at multiple levels
     *
     * @param added non-null change listener
     */
    @Override
    public final void addSpectrumHolderListener(SpectrumHolderListener added) {
        if (!spectrumHolderListeners.contains(added))
            spectrumHolderListeners.add(added);
        SpectrumUtilities.guaranteeSerializable(this);   // todo remove
    }

    /**
     * remove a change listener
     *
     * @param removed non-null change listener
     */
    @Override
    public final void removeSpectrumHolderListener(SpectrumHolderListener removed) {
        while (spectrumHolderListeners.contains(removed))
            spectrumHolderListeners.remove(removed);
    }


    /**
     * notify any state change listeners - probably should
     * be protected but is in the interface to form an event cluster
     */
    protected void notifySpectrumHolderListeners(boolean isAdd, ISpectrum... spectra) {
        if (spectrumHolderListeners.isEmpty())
            return;
        for (SpectrumHolderListener listener : spectrumHolderListeners) {
            if (isAdd)
                listener.onSpectraAdd(this, spectra);
            else
                listener.onSpectraRemove(this, spectra);
        }
    }

    /**
     * if possible use the highest
     *
     * @return The cluster's id
     */
    @Override
    public String getId() {
        if (id == null) {
            // if no id is set, generate a new unique id
            id = UUID.randomUUID().toString();
        }
        return id;
    }

    @Override
    public String getSpectralId() {
        StringBuilder sb = new StringBuilder();
        List<String> spectralIds = new ArrayList<String>(getSpectralIds());

        if (spectralIds.size() > 1) {
            Collections.sort(spectralIds);
            sb.append(spectralIds.get(0));
            for (int i = 1; i < spectralIds.size(); i++) {
                sb.append(",");
                sb.append(spectralIds.get(i));
            }
            return sb.toString();
        } else {
            return spectralIds.get(0);
        }
    }

    public void setId(String id) {
        this.id = id;
    }

    @Override
    public float getPrecursorMz() {
        ISpectrum consensusSpectrum1 = getConsensusSpectrum();
        if (consensusSpectrum1 == null)
            return 0;
        return consensusSpectrum1.getPrecursorMz();
    }

    @Override
    public int getPrecursorCharge() {
        ISpectrum consensusSpectrum1 = getConsensusSpectrum();
        if (consensusSpectrum1 == null)
            return 0;
        return consensusSpectrum1.getPrecursorCharge();
    }

    /**
     * all internally spectrum
     */
    @Override
    public List<ISpectrum> getHighestQualitySpectra() {
        throw new UnsupportedOperationException("Marked deprecated.");
    }


    /**
     * real spectrum with the highest quality - this is a
     * good way to compare clusters
     *
     * @return !null spectrum
     */
    @Override
    public ISpectrum getHighestQualitySpectrum() {
        return getConsensusSpectrum();
    }


    @Override
    public ISpectrum getConsensusSpectrum() {
        return consensusSpectrumBuilder.getConsensusSpectrum();
    }

    @Override
    public List<ISpectrum> getClusteredSpectra() {
        return new ArrayList<ISpectrum>(clusteredSpectra);
    }

    @Override
    public IConsensusSpectrumBuilder getConsensusSpectrumBuilder() {
        return consensusSpectrumBuilder;
    }

    @Override
    public int getClusteredSpectraCount() {
        return clusteredSpectra.size();
    }

    @Override
    public void addSpectra(ISpectrum... merged) {
        if (merged != null && merged.length > 0) {
            boolean spectrumAdded = false;
            final ArrayList<ISpectrum> added = new ArrayList<ISpectrum>();

            for (ISpectrum spectrumToMerge : merged) {
                // ignore spectra that have already been added
                if (spectraIds.contains(spectrumToMerge.getId()))
                    continue;

                spectraIds.add(spectrumToMerge.getId());

                spectrumAdded = true;
                ISpectrum spectrumWithoutPeaks = new Spectrum(spectrumToMerge, Collections.EMPTY_LIST);
                // only store spectra without peak lists to save memory
                clusteredSpectra.add(spectrumWithoutPeaks);
                added.add(spectrumToMerge);
            }

            if (spectrumAdded)
                notifySpectrumHolderListeners(true, added.toArray(new ISpectrum[added.size()]));   // tell other interested parties  true says this is an add
        }
    }

    /**
     * This function enables the merging of clusters that do not save peak lists
     *
     * @param cluster An ICluster to add to the cluster
     */
    public void addCluster(ICluster cluster) {
        if (cluster.storesPeakLists()) {
            ISpectrum[] spectraToAdd = new ISpectrum[cluster.getClusteredSpectraCount()];
            addSpectra(cluster.getClusteredSpectra().toArray(spectraToAdd));
        } else {
            // it doesn't store peak lists so only merge the consensus spectra
            if (cluster.getClusteredSpectra() == null || cluster.getClusteredSpectra().isEmpty())
                return;

            // simply add the consensus spectrum to this one
            consensusSpectrumBuilder.addConsensusSpectrum(cluster.getConsensusSpectrumBuilder());

            // add the spectra and their ids
            clusteredSpectra.addAll(cluster.getClusteredSpectra());

            // save the spectra ids
            for (ISpectrum spectrum : cluster.getClusteredSpectra())
                spectraIds.add(spectrum.getId());

            notifySpectrumHolderListeners(true, cluster.getClusteredSpectra().toArray(new ISpectrum[cluster.getClusteredSpectra().size()]));   // tell other interested parties  true says this is an add
        }
    }


    /**
     * stable clusters do not support remove others do
     *
     * @return as above
     */
    @Override
    public boolean isRemoveSupported() {
        return false;
    }

    @Override
    public void removeSpectra(ISpectrum... removed) {
        throw new UnsupportedOperationException("Remove not supported");
    }

    /**
     * return a property of null if none exists
     * See ISpectrum for known property names
     *
     * @param key String representing the property's name
     * @return possible null value
     */
    @Override
    public String getProperty(String key) {
        return properties.getProperty(key);
    }


    /**
     * @param key
     * @param value
     */
    @Override
    public void setProperty(String key, String value) {
        if (key == null)
            return;
        if (value == null) {
            properties.remove(key);
            return;
        }

        properties.setProperty(key, value);
    }

    /**
     * Only for internal use in copy constructor
     * Note this is not safe
     * This is not really deprecated but it warns only for
     * internal use
     */
    @Override
    public Properties getProperties() {
        return properties;
    }

    /**
     * sort by mz - might be useful
     * NOTE we cannot use ClusterComparator because we want to compare spectralIds
     *
     * @param o The cluster to compare to
     * @return Integer representing the sort order
     */
    @Override
    public int compareTo(ICluster o) {
        if (o == this)
            return 0;
        try {
            int ret = CompareTo.compare(getPrecursorMz(), o.getPrecursorMz());
            if (ret != 0)
                return ret;

            if (o.getClusteredSpectraCount() != getClusteredSpectraCount()) {
                return getClusteredSpectraCount() < o.getClusteredSpectraCount() ? -1 : 1;
            }

            String spectra = getSpectralId();
            String otherSpectra = o.getSpectralId();
            return spectra.compareTo(otherSpectra);

        } catch (IllegalStateException e) {
            //  give up use hash code
        }

        int hash1 = hashCode();
        int hash2 = o.hashCode();
        if (hash1 != hash2)
            return hash1 < hash2 ? -1 : 1;

        return 0;
    }


    /**
     * like equals but weaker - says other is equivalent to this
     *
     * @param o poiibly null other object
     * @return true if other is "similar enough to this"
     */
    @Override
    public boolean equivalent(ICluster o) {
        if (o == this)
            return true;
        if (getPrecursorCharge() != o.getPrecursorCharge())
            return false;
        double del = o.getPrecursorMz() - getPrecursorMz();
        double abs = Math.abs(del);
        if (abs > MZIntensityUtilities.SMALL_MZ_DIFFERENCE) {
            return false;
        }

        List<ISpectrum> spc1 = clusteredSpectra;
        List<ISpectrum> spc2 = o.getClusteredSpectra();

        if (spc1.size() != spc2.size()) {
            return false;
        }
        // NOTE for only one spectrum in the cluster you are really comparing the only spectrum
        if (spc1.size() <= 1) {

            final ISpectrum consensusSpectrum = getConsensusSpectrum();
            final ISpectrum oConsensusSpectrum = o.getConsensusSpectrum();

            return consensusSpectrum.equivalent(oConsensusSpectrum);

        } else {
            for (int i = 0; i < spc1.size(); i++) {
                ISpectrum pk1 = spc1.get(i);
                ISpectrum pk2 = spc2.get(i);
                boolean test = !pk2.equivalent(pk1);
                if (test)
                    return false;
            }
            return true;
        }
    }

    @Override
    public String toString() {
        double precursorMZ = getPrecursorMz();
        String text =
                "charge= " + getPrecursorCharge() + "," +
                        "mz= " + String.format("%10.3f", precursorMZ).trim() + "," +
                        "count= " + clusteredSpectra.size() +
                        ", spectrum = ";
        for (ISpectrum s : clusteredSpectra)
            text += s.getId() + ",";

        text = text.substring(0, text.length() - 1);
        return text;
    }

    /**
     * Saves the comparison match in the best matches array
     *
     * @param id Id of the cluster that the comparison was performed with
     * @param similarity The similarity score to store for this comparison
     */
    public void saveComparisonResult(String id, float similarity) {
        if (bestComparisonMatches.size() >= SAVED_COMPARISON_MATCHES && similarity < lowestBestComparisonSimilarity)
            return;

        ComparisonMatch comparisonMatch = new ComparisonMatch(id, similarity);
        bestComparisonMatches.add(comparisonMatch);

        // remove the lowest matches if necessary
        if (bestComparisonMatches.size() >= SAVED_COMPARISON_MATCHES) {
            Collections.sort(bestComparisonMatches);
            int tooManyMatches = bestComparisonMatches.size() - SAVED_COMPARISON_MATCHES;
            for (int i = 0; i < tooManyMatches; i++)
                bestComparisonMatches.remove(0); // natural order is lowest to highest, remove lowest

            lowestBestComparisonSimilarity = bestComparisonMatches.get(0).getSimilarity();
        }
        // if nothing is removed check whether the new lowest similarity is the lowest
        else if (similarity < lowestBestComparisonSimilarity) {
            lowestBestComparisonSimilarity = similarity;

        }

        bestComparisonMatchIds = null; // delete to mark as dirty
    }

    /**
     * Checks whether a given spectrum id is part of the
     * best similarity matches.
     *
     * @param id The other cluster's id
     * @return Boolean indicating whether the comparison scored among the top N
     */
    public boolean isInBestComparisonResults(String id) {
        if (bestComparisonMatchIds == null) {
            bestComparisonMatchIds = new HashSet<String>(bestComparisonMatches.size());

            for (ComparisonMatch comparisonMatch : bestComparisonMatches)
                bestComparisonMatchIds.add(comparisonMatch.getSpectrumId());
        }

        return bestComparisonMatchIds.contains(id);
    }

    @Override
    public boolean storesPeakLists() {
        return false;
    }

    @Override
    public List<ComparisonMatch> getComparisonMatches() {
        return Collections.unmodifiableList(bestComparisonMatches);
    }

    @Override
    public void setComparisonMatches(List<ComparisonMatch> comparisonMatches) {
        this.bestComparisonMatches.clear();
        if (comparisonMatches != null && comparisonMatches.size() > 0) {
            this.bestComparisonMatches.addAll(comparisonMatches);

            Collections.sort(bestComparisonMatches);
            lowestBestComparisonSimilarity = bestComparisonMatches.get(0).getSimilarity();
        } else {
            lowestBestComparisonSimilarity = 0;
        }

        bestComparisonMatchIds = null; // delete to mark as dirty
    }

    @Override
    public boolean isKnownComparisonMatch(String clusterId) {
        if (bestComparisonMatches.size() == 0)
            return false;

        for (ComparisonMatch comparisonMatch : bestComparisonMatches) {
            if (comparisonMatch.getSpectrumId().equals(clusterId))
                return true;
        }

        return false;
    }
}

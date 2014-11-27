package uk.ac.ebi.pride.spectracluster.cluster;

import uk.ac.ebi.pride.spectracluster.consensus.IConsensusSpectrumBuilder;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.util.*;

import java.util.*;
import java.util.concurrent.CopyOnWriteArrayList;


/**
 * Default implementation of ISpectralCluster
 *
 * @author Steve Lewis
 * @author Rui Wang
 * @version $Id$
 */
public class SpectralCluster implements ICluster {

    private String id;
    // holds a list of the top  SpectralQualityHolder.NUMBER_SPECTRA_FOR_CONSENSUS = 20;
    // quality spectra - these can be use to build a concensus of quality
    // Note all adds and removes are done by registering as a SpectrumHolderListener
    private final SpectralQualityHolder qualityHolder;
    private final List<SpectrumHolderListener> spectrumHolderListeners = new CopyOnWriteArrayList<SpectrumHolderListener>();
    private final Set<String> spectraIds = new HashSet<String>();
    private final Properties properties = new Properties();

    private final List<ISpectrum> clusteredSpectra = new ArrayList<ISpectrum>();
    private final IConsensusSpectrumBuilder consensusSpectrumBuilder;

    public SpectralCluster(ICluster copied, IConsensusSpectrumBuilder consensusSpectrumBuilder) {
        this(copied.getId(), consensusSpectrumBuilder);

        final List<ISpectrum> clusteredSpectra1 = copied.getClusteredSpectra();
        addSpectra(clusteredSpectra1.toArray(new ISpectrum[clusteredSpectra1.size()]));
    }


    public SpectralCluster(String id, IConsensusSpectrumBuilder consensusSpectrumBuilder) {
        this.id = id;
        this.consensusSpectrumBuilder = consensusSpectrumBuilder;
        addSpectrumHolderListener(this.consensusSpectrumBuilder);
        this.qualityHolder = new SpectralQualityHolder();
        addSpectrumHolderListener(qualityHolder);
        SpectrumUtilities.guaranteeSerializable(this);   // todo remove
    }

    /**
     * return a set of all ids
     *
     * @return
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
     * @return
     */
    @Override
    public String getId() {
        // in unstable clusters use id of the highest quality spectrum
        if (id == null) {
            id = getSpectralId();
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
        return qualityHolder.getHighestQualitySpectra();
    }


    /**
     * real spectrum with the highest quality - this is a
     * good way to compare clusters
     *
     * @return !null spectrum
     */
    @Override
    public ISpectrum getHighestQualitySpectrum() {
        return qualityHolder.getHighestQualitySpectrum();
    }


    @Override
    public ISpectrum getConsensusSpectrum() {
        //  todo I think in a cluster with one spectrum we do not need to build a
        //   consensus spectrum
        if(clusteredSpectra.size() == 1)
            return clusteredSpectra.get(0);
        return consensusSpectrumBuilder.getConsensusSpectrum();
    }

    @Override
    public List<ISpectrum> getClusteredSpectra() {
        return new ArrayList<ISpectrum>(clusteredSpectra);
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
                spectraIds.add(spectrumToMerge.getId());
                if (!clusteredSpectra.contains(spectrumToMerge)) {
                    spectrumAdded = true;
                    clusteredSpectra.add(spectrumToMerge);
                    added.add(spectrumToMerge);
                }
            }
            if (spectrumAdded)
                notifySpectrumHolderListeners(true, added.toArray(new ISpectrum[added.size()]));   // tell other interested parties  true says this is an add
        }
        SpectrumUtilities.guaranteeSerializable(this);   // todo remove
    }


    /**
     * stable clusters do not support remove others do
     *
     * @return as above
     */
    @Override
    public boolean isRemoveSupported() {
        // todo Johannes stable clusters do not support removal - do we????
        // return !isStable();
        return true;
    }

    @Override
    public void removeSpectra(ISpectrum... removed) {
        if (!isRemoveSupported())
            throw new UnsupportedOperationException("Remove not supported");

        if (removed != null && removed.length > 0) {
            for (ISpectrum spectrumToRemove : removed) {
                spectraIds.remove(spectrumToRemove.getId());
                clusteredSpectra.remove(spectrumToRemove);
            }

            notifySpectrumHolderListeners(false, removed); // tell other interested parties  false says this is a remove
        }
    }

    /**
     * return a property of null if none exists
     * See ISpectrum for known property names
     *
     * @param key
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
        if(key == null)
            return;
        if( value == null)   {
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
     *   NOTE we cannot use ClusterComparator because we want to compare spectralIds
     * @param o
     * @return
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
}

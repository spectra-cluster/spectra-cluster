package uk.ac.ebi.pride.spectracluster.consensus;

import uk.ac.ebi.pride.spectracluster.cluster.ISpectrumHolder;
import uk.ac.ebi.pride.spectracluster.cluster.SpectrumHolderListener;
import uk.ac.ebi.pride.spectracluster.filter.BinnedHighestNPeakFilter;
import uk.ac.ebi.pride.spectracluster.filter.IPeakFilter;
import uk.ac.ebi.pride.spectracluster.normalizer.TotalIntensityNormalizer;
import uk.ac.ebi.pride.spectracluster.spectrum.IPeak;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.spectrum.Peak;
import uk.ac.ebi.pride.spectracluster.spectrum.Spectrum;
import uk.ac.ebi.pride.spectracluster.util.Defaults;
import uk.ac.ebi.pride.spectracluster.util.MZIntensityUtilities;
import uk.ac.ebi.pride.spectracluster.util.PeakUtilities;
import uk.ac.ebi.pride.spectracluster.util.comparator.PeakIntensityComparator;
import uk.ac.ebi.pride.spectracluster.util.comparator.PeakMzComparator;

import java.util.*;

/**
 * This is a rewrite of the original FrankEtAlConsensuSpectrumBuilder and produces nearly identical results
 * as the original class. The differences that are currently observed only occur under special circumstances and
 * are caused by rounding differences between double and float. This may lead to different peak merging results.
 * <p/>
 * The main difference of this class to the original one is that the consensus spectrum is generated from one
 * internal "crowded" spectrum. This "crowded" spectrum contains all peaks of all added spectra. Thereby, when new
 * spectra are added or removed it is not necessary to reprocess all previous spectra.
 * <p/>
 * Theoretically, memory can be saved since identical peaks between different spectra are only stored once. This effect
 * can be increased by rounding the m/z value to f.e. 3 digits after the comma (this does not influence the final
 * consensus spectrum). Since this change increased the processing time of adding spectra at least two fold, it was
 * removed again. The function "round" that was used to the rounding is still part of the class.
 * NOTE This is Johannes code refactored so that critical methods take and return peak lists rather than
 * using an internal list - this allows steps to be individually rewritten and tested
 * Also internal methods are protected allowing tests to access them
 *
 * todo: development
 *
 */
@Deprecated
public class FilteredConsensusSpectrum implements IConsensusSpectrumBuilder {

    public static final ConcensusSpectrumBuilderFactory FACTORY = new ConsensusSpectrumFactory();

    /**
     * always use the factory to get an instance
     */
    public static class ConsensusSpectrumFactory implements ConcensusSpectrumBuilderFactory {
        private ConsensusSpectrumFactory() {
        }

        /**
         * build a new instance of the cpectrum builder
         *
         * @return !null instance
         */
        @Override
        public IConsensusSpectrumBuilder getConsensusSpectrumBuilder() {
            return new FilteredConsensusSpectrum();
        }
    }


    public static final int DEFAULT_PEAKS_TO_KEEP = 5;
    public static final int SIZE_TO_ADD_EVERY_TIME = 100;   // if less then this add all peaks
    public static final float FRACTION_OF_LOWEST_PEAK_TOKEEP = 0.40F; // do not keep peaks this much smaller than what er currently keep


    public static final float NOISE_FILTER_INCREMENT = 100;

    private final IPeakFilter filter;
    private final String id;
    protected int nSpectra;
    protected boolean isDirty;
    protected float averagePrecursorMz;
    protected float sumPrecursorMz;
    protected float averagePrecursorIntens;
    protected float sumPrecursorIntens;
    protected float lowestConcensusPeak;
    protected int averageCharge;
    protected int sumCharge;
    protected ISpectrum consensusSpectrum;
    protected final List<SpectrumHolderListener> listeners = new ArrayList<SpectrumHolderListener>();

    protected final String methodName = "Crowded Consensus Spectrum Builder";
    protected final String methodVersion = "0.1";

    /**
     * The m/z threshold to consider two peaks identical
     */
    protected static final float FINAL_MZ_THRESHOLD = 0.4F;
    /**
     * The m/z threshold for identical peaks is not applied instantly, but gradually increased
     * to reach the final threshold. Each iteration increases the threshold by MZ_THRESHOLD_STEP.
     */
    protected static final float MZ_THRESHOLD_STEP = 0.1F;
    /**
     * Defines whether m/z values should be rounded. Thereby, more peaks are considered identical
     * without reducing accuracy.
     */
    public final static boolean USE_ROUNDING = false;
    /**
     * Defines the precision used if rounding is enabled.
     */
    public final static int MZ_PRECISION = 1000;

    /**
     * Holds all peaks from all added spectra. In case an exact m/z is found twice, the intensities are added.
     * The array must always be sorted according to m/z.
     */
    private final List<IPeak> allPeaks = new ArrayList<IPeak>();
    /**
     * Holds all peaks that are waiting to be added but have not been
     */
    private final Set<IPeak> heldPeaks = new HashSet<IPeak>();
    /**
     * Peaks of the actual consensusSpectrum
     */
    private final List<IPeak> consensusPeaks = new ArrayList<IPeak>();


    /**
     * private to force use of the factory
     */
    private FilteredConsensusSpectrum() {
        this(null);
    }

    /**
     * private to force use of the factory
     */
    private FilteredConsensusSpectrum(String id) {

        this.id = id;
        filter = BinnedHighestNPeakFilter.DEFAULT;
    }


    /**
     * expose a normally private field for testing
     *
     * @return
     */
    @SuppressWarnings("UnusedDeclaration")
    protected List<IPeak> getInternalPeaks() {
        return consensusPeaks;
    }


    @Override  // TODO JG this class only correctly supports normalized spectra. Make sure the spectra are normalized
    public void addSpectra(ISpectrum... merged) {
        if (merged.length < 1)
            return;

        ISpectrum[] filtered = new ISpectrum[merged.length];

        for (int i = 0; i < merged.length; i++) {
            ISpectrum spectrum = merged[i];

            List<IPeak> oldSpectrumPeaks = spectrum.getPeaks();

            // filter the peaks by bin
            List<IPeak> filteredPeaks = filter.filter(oldSpectrumPeaks);
            filteredPeaks = TotalIntensityNormalizer.DEFAULT.normalizePeaks(filteredPeaks);
            // add the filtered peaks
            addPeaks(filteredPeaks);

            filtered[i] = new Spectrum(spectrum, filteredPeaks);
            sumCharge += spectrum.getPrecursorCharge();
            sumPrecursorMz += spectrum.getPrecursorMz();
            sumPrecursorIntens += 0;

            nSpectra++;


        }

        setIsDirty(true);

        for (SpectrumHolderListener listener : listeners)
            listener.onSpectraAdd(this, filtered);
    }

    @Override
    public void removeSpectra(ISpectrum... removed) {
        if (removed.length < 1)
            return;

        for (ISpectrum spectrum : removed) {
            List<IPeak> spectrumPeaks = spectrum.getPeaks();
            removePeaks(spectrumPeaks);

            sumCharge -= spectrum.getPrecursorCharge();
            sumPrecursorMz -= spectrum.getPrecursorMz();
            sumPrecursorIntens -= 0;

            nSpectra--;
        }

        setIsDirty(true);

        for (SpectrumHolderListener listener : listeners)
            listener.onSpectraRemove(this, removed);
    }

    /**
     * stable clusters do not support remove others do
     *
     * @return as above
     */
    @Override
    public boolean isRemoveSupported() {
        return true;
    }

    /**
     * Removes the passed peaks from the "crowded" spectrum
     * allPeaks.
     *
     * @param peaksToRemove
     */
    protected void removePeaks(List<IPeak> peaksToRemove) {
        int posAllPeaks = 0;
        //noinspection ForLoopReplaceableByForEach
        for (int i = 0; i < peaksToRemove.size(); i++) {
            IPeak peakToRemove = peaksToRemove.get(i);
            double mzToRemove = peakToRemove.getMz();

            if (USE_ROUNDING)
                mzToRemove = MZIntensityUtilities.round(mzToRemove, MZ_PRECISION);

            for (int j = posAllPeaks; j < allPeaks.size(); j++) {
                IPeak currentExistingPeak = allPeaks.get(j);

                if (mzToRemove < currentExistingPeak.getMz()) {
                    // TODO @Rui/Steve: This means that the peak does not exist, should we throw an exception here?
                    posAllPeaks = j;
                    break;
                }

                if (mzToRemove == currentExistingPeak.getMz()) {
                    allPeaks.set(j, new Peak(
                                    currentExistingPeak.getMz(),
                                    currentExistingPeak.getIntensity() - peakToRemove.getIntensity(),
                                    currentExistingPeak.getCount() - 1)
                    );

                    posAllPeaks = j;
                    break;
                }
            }
        }

        // clear all peaks with count < 1
        List<IPeak> tmp = new ArrayList<IPeak>();
        for (IPeak p : allPeaks) {
            if (p.getCount() > 0)
                tmp.add(p);
        }
        allPeaks.clear();
        allPeaks.addAll(tmp);
    }

    @Override
    public void addSpectrumHolderListener(SpectrumHolderListener added) {
        listeners.add(added);
    }

    @Override
    public void removeSpectrumHolderListener(SpectrumHolderListener removed) {
        listeners.remove(removed);
    }

    /**
     * Modified to slow peek adding as the number of spectra gets large
     *
     * @param peaksToAdd
     */
    protected void addPeaks(List<IPeak> peaksToAdd) {
        if (nSpectra < SIZE_TO_ADD_EVERY_TIME) {
            internalAddPeaks(peaksToAdd);
            return;
        }
        storeHeldPeaks(peaksToAdd);
        if (nSpectra < 10 * SIZE_TO_ADD_EVERY_TIME) {
            if (nSpectra % SIZE_TO_ADD_EVERY_TIME == 0) {
                addHeldPeaks();
            }
            return;
        }
        if (nSpectra % 4 * SIZE_TO_ADD_EVERY_TIME == 0) {
            addHeldPeaks();
        }

    }


    /**
     * add all peaks being held
     */
    protected void storeHeldPeaks(List<IPeak> peaksToAdd) {
        if (nSpectra < SIZE_TO_ADD_EVERY_TIME)
            throw new IllegalStateException("cannot add without a fairly large cluster");
        // force one computation if we have never done this
        if (lowestConcensusPeak == 0)
            getConsensusSpectrum();

        //noinspection UnnecessaryLocalVariable,UnusedDeclaration,UnusedAssignment
        int skipped = 0;
        float minimumKeptPeak = lowestConcensusPeak * FRACTION_OF_LOWEST_PEAK_TOKEEP;
        for (IPeak iPeak : peaksToAdd) {
            if (iPeak.getIntensity() > minimumKeptPeak) {
                heldPeaks.add(iPeak);
            } else {
                skipped++;
            }
        }

    }


    /**
     * add all peaks being held
     */
    protected void addHeldPeaks() {
        List<IPeak> added = new ArrayList<IPeak>(heldPeaks);
        Collections.sort(added);
        internalAddPeaks(added);
        heldPeaks.clear();
    }


    /**
     * Adds the passed peaks to the "crowded" internal spectrum (allPeaks). The precursor m/z
     * values are rounded to MZ_PRECISION digist after the comma. This increases the probability that
     * two peaks have the identical precursor m/z and only have to be stored as one peak.
     *
     * @param peaksToAdd
     */
    protected void internalAddPeaks(List<IPeak> peaksToAdd) {
        //TODO @jg: build in a check to find if peaks are not sorted according to m/z
        int posAllPeaks = 0;
        List<IPeak> newPeaks = new ArrayList<IPeak>(); // peaks with m/z values that do not yet exist

        //noinspection ForLoopReplaceableByForEach
        for (int i = 0; i < peaksToAdd.size(); i++) {
            IPeak peakToAdd = peaksToAdd.get(i);
            double mzToAdd = peakToAdd.getMz();

            if (USE_ROUNDING)
                mzToAdd = MZIntensityUtilities.round(mzToAdd, MZ_PRECISION);

            boolean wasAdded = false;

            for (int j = posAllPeaks; j < allPeaks.size(); j++) {
                IPeak currentExistingPeak = allPeaks.get(j);

                if (mzToAdd < currentExistingPeak.getMz()) {
                    newPeaks.add(new Peak((float) mzToAdd, peakToAdd.getIntensity(), peakToAdd.getCount()));
                    posAllPeaks = j;
                    wasAdded = true;
                    break;
                }

                if (mzToAdd == currentExistingPeak.getMz()) {
                    allPeaks.set(j, new Peak(
                                    currentExistingPeak.getMz(),
                                    peakToAdd.getIntensity() + currentExistingPeak.getIntensity(),
                                    currentExistingPeak.getCount() + peakToAdd.getCount())
                    );
                    posAllPeaks = j;
                    wasAdded = true;
                    break;
                }
            }

            if (!wasAdded)
                newPeaks.add(new Peak((float) mzToAdd, peakToAdd.getIntensity(), peakToAdd.getCount()));
        }

        // add all new peaks
        allPeaks.addAll(newPeaks);
        Collections.sort(allPeaks, new PeakMzComparator());
    }

    /**
     * Indicates whether the current data is up-to-date. As soon as spectra are
     * removed or added the spectrum is considered dirty.
     *
     * @return
     */
    protected boolean isDirty() {
        return isDirty;
    }

    protected void setIsDirty(boolean isDirty) {
        this.isDirty = isDirty;
    }

    /**
     * Updates all properties of the consensus spectrum as well as the actual consensus
     * spectrum.
     */
    protected void update() {
        if (nSpectra > 0) {
            averagePrecursorMz = sumPrecursorMz / nSpectra;
            averageCharge = sumCharge / nSpectra;
            averagePrecursorIntens = sumPrecursorIntens / nSpectra;
        } else {
            averagePrecursorMz = 0;
            averageCharge = 0;
            averagePrecursorIntens = 0;
        }


        if (allPeaks.size() < 1) {
            List<IPeak> empty = new ArrayList<IPeak>();
            consensusSpectrum = new Spectrum(id, averageCharge, averagePrecursorMz, Defaults.getDefaultQualityScorer(), Collections.EMPTY_LIST);
            setIsDirty(false);
            return;
        }
        List<IPeak> newPeaks = findConsensusPeaks(allPeaks, DEFAULT_PEAKS_TO_KEEP, nSpectra);

        // update the consensus spectrum
        consensusPeaks.clear();
        // allPeaks is always sorted according to precursor m/z
        consensusPeaks.addAll(newPeaks);

        // find the minimum peak - much smaller and we dont need to save
        float minimumConsensusPeak = Float.MAX_VALUE;
        for (IPeak peak : consensusPeaks) {
            final float intensity = peak.getIntensity();
            if (intensity < minimumConsensusPeak && intensity > 0)
                minimumConsensusPeak = intensity;
        }
        lowestConcensusPeak = minimumConsensusPeak;

        // create the ConsensusSpectrum object
        consensusSpectrum = new Spectrum(id, averageCharge, averagePrecursorMz, Defaults.getDefaultQualityScorer(), consensusPeaks);

        setIsDirty(false);
    }

    /**
     * refactored to make testing easier
     *
     * @param input !null set of all peaks
     * @return !null set of  consensus peaks
     */
    protected static List<IPeak> findConsensusPeaks(List<IPeak> input, int peaksToKeep, int nSpectra) {


        // Step 1: merge identical peaks
        List<IPeak> ret = mergeIdenticalPeaks(input);


        // Step 2: adapt the peak intensities based on the probability that the peak has been observed
        ret = adaptPeakIntensities(ret, nSpectra);

        // Step 3: filter the spectrum
        ret = filterNoise(ret, peaksToKeep);
        return ret;
    }

    /**
     * Filters the consensus spectrum keeping only the top 5 peaks per 100 m/z
     */
    protected static List<IPeak> filterNoise(List<IPeak> inp, int peaksInBinToKeep) {
        List<IPeak> filteredSpectrum = new ArrayList<IPeak>();

        int lowerBound = 0;
        // process the peaks using a sliding window of 100 m/z
        for (double startMz = 0, endMz = NOISE_FILTER_INCREMENT; endMz <= MZIntensityUtilities.HIGHEST_USABLE_MZ; endMz += NOISE_FILTER_INCREMENT, startMz += NOISE_FILTER_INCREMENT) {
            List<IPeak> peakBuffer = new ArrayList<IPeak>();

            // set the lower bound
            for (int i = lowerBound; i < inp.size(); i++) {
                if (inp.get(i).getMz() >= startMz) {
                    lowerBound = i;
                    break;
                }
            }

            if (inp.get(lowerBound).getMz() < startMz)
                continue;

            for (int i = lowerBound; i < inp.size(); i++) {
                if (inp.get(i).getMz() <= endMz) {
                    peakBuffer.add(inp.get(i));
                } else {
                    lowerBound = i;
                    break;
                }
            }

            if (peakBuffer.size() < 1)
                continue;

            Collections.sort(peakBuffer, PeakIntensityComparator.INSTANCE);

            List<IPeak> fivePeaks = new ArrayList<IPeak>(peaksInBinToKeep);

            for (int i = 0; i < peaksInBinToKeep && i < peakBuffer.size(); i++)
                fivePeaks.add(peakBuffer.get(i));

            Collections.sort(fivePeaks, new PeakMzComparator());
            filteredSpectrum.addAll(fivePeaks);
        }

        return filteredSpectrum;
    }

    /**
     * Adapt the peak intensities in consensusPeaks using the following formula:
     * I = I * (0.95 + 0.05 * (1 + pi))^5
     * where pi is the peaks probability
     */
    protected static List<IPeak> adaptPeakIntensities(List<IPeak> inp, int nSpectra) {

        int originalCount = PeakUtilities.getTotalCount(inp);   // for debugging

        List<IPeak> ret = new ArrayList<IPeak>(inp);
        for (int i = 0; i < ret.size(); i++) {
            IPeak peak = ret.get(i);
            float peakProbability = (float) peak.getCount() / (float) nSpectra;
            float newIntensity = (float) (peak.getIntensity() * (0.95 + 0.05 * Math.pow(1 + peakProbability, 5)));

            ret.set(i, new Peak(peak.getMz(), newIntensity, peak.getCount()));
        }

        int finalCount = PeakUtilities.getTotalCount(ret);   // for debugging
        if (originalCount != finalCount)
            throw new IllegalStateException("Peak merge changed total count");

        return ret;
    }

    /**
     * Merges identical peaks in the consensusPeaks List based on FINAL_MZ_THRESHOLD and
     * MZ_THRESHOLD_STEP.
     */
    protected static List<IPeak> mergeIdenticalPeaks(List<IPeak> inPeaks) {
        int originalCount = PeakUtilities.getTotalCount(inPeaks);   // for debugging

        List<IPeak> ret = new ArrayList<IPeak>();
        if (inPeaks.size() == 0)
            return ret; // should never happen

        for (float range = MZ_THRESHOLD_STEP; range <= FINAL_MZ_THRESHOLD; range += MZ_THRESHOLD_STEP) {
            List<IPeak> newPeakList = new ArrayList<IPeak>();
            IPeak currentPeak = inPeaks.get(0);

            for (int i = 1; i < inPeaks.size(); i++) {
                IPeak nextPeak = inPeaks.get(i);

                // check whether the next peak should be considered identical to the current one
                final float nextPeakMz = nextPeak.getMz();
                final float currentPeakMz = currentPeak.getMz();
                final float testLimit = currentPeakMz + range;
                if (nextPeakMz <= testLimit) {
                    // calculate the new weighted m/z
                    final double nextPeakIntensity = nextPeak.getIntensity();
                    final double currentPeakIntensity = currentPeak.getIntensity();
                    final double totalIntensity = nextPeakIntensity + currentPeakIntensity;
                    final double nextPeakFraction = nextPeakIntensity / totalIntensity;
                    final double currentPeakFraction = currentPeakIntensity / totalIntensity;

                    //noinspection UnnecessaryLocalVariable,UnusedDeclaration,UnusedAssignment
                    final double totalFraction = nextPeakFraction + currentPeakFraction;


                    double weightedMz = (nextPeakFraction * nextPeakMz) + (currentPeakFraction * currentPeakMz);


                    final double intensity = currentPeakIntensity + nextPeakIntensity;
                    final int count = currentPeak.getCount() + nextPeak.getCount();
                    //noinspection UnnecessaryLocalVariable,UnusedDeclaration,UnusedAssignment
                    IPeak newPeak = new Peak((float) weightedMz, (float) intensity, count);

                    currentPeak = newPeak;
                } else {
                    // by adding the peak in the else clause, peaks that were merged are not included in the new Peak
                    // list and are thereby removed from the consensusPeaks
                    newPeakList.add(currentPeak);
                    currentPeak = nextPeak;
                }
            }
            newPeakList.add(currentPeak);

            ret.clear();
            ret.addAll(newPeakList);

        }

        int finalCount = PeakUtilities.getTotalCount(ret);   // for debugging
        if (originalCount != finalCount)
            throw new IllegalStateException("Peak merge changed total count");

        return ret;
    }

    @Override
    public ISpectrum getConsensusSpectrum() {
        if (isDirty())
            update();

        return internalGetConcensusSpectrum();
    }

    /**
     * access for test purposes
     *
     * @return possibly null spectrum
     */
    protected ISpectrum internalGetConcensusSpectrum() {
        return consensusSpectrum;
    }


    @Override
    public void clear() {
        sumCharge = 0;
        sumPrecursorMz = 0;
        sumPrecursorIntens = 0;
        nSpectra = 0;

        allPeaks.clear();
        consensusPeaks.clear();

        setIsDirty(true);
    }

    @Override
    public int getSpectraCount() {
        return nSpectra;
    }

    @Override
    public String getName() {
        return methodName;
    }

    @Override
    public String getCurrentVersion() {
        return methodVersion;
    }

    @Override
    public void onSpectraAdd(ISpectrumHolder holder, ISpectrum... added) {
        addSpectra(added);
    }

    @Override
    public void onSpectraRemove(ISpectrumHolder holder, ISpectrum... removed) {
        removeSpectra(removed);
    }
}

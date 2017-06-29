package uk.ac.ebi.pride.spectracluster.consensus;

import uk.ac.ebi.pride.spectracluster.cluster.ISpectrumHolder;
import uk.ac.ebi.pride.spectracluster.cluster.SpectrumHolderListener;
import uk.ac.ebi.pride.spectracluster.spectrum.IPeak;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.spectrum.Peak;
import uk.ac.ebi.pride.spectracluster.spectrum.Spectrum;
import uk.ac.ebi.pride.spectracluster.util.Defaults;
import uk.ac.ebi.pride.spectracluster.util.comparator.PeakMzComparator;
import uk.ac.ebi.pride.spectracluster.util.function.IFunction;
import uk.ac.ebi.pride.spectracluster.util.function.peak.BinnedHighestNPeakFunction;

import java.util.*;

/**
 * This is a greedy version of the FrankEtAlConsensusSpectrumBuilder. It only
 * supports the addition of spectra but not their removal. Thereby, the original
 * peaks do not have to be kept.
 *
 * @author Johannes Griss
 */
public class GreedyConsensusSpectrum implements IConsensusSpectrumBuilder {
    /**
     * Peaks to keep per 100 m/z during noise filtering
     */
    public static final int DEFAULT_PEAKS_TO_KEEP = 5;
    /**
     * Sliding window size to use when applying the noise filter
     */
    public static final float NOISE_FILTER_INCREMENT = 100;
    /**
     * The filter to use for (noise) filtering of the final consensus spectrum
     */
    private final static IFunction<List<IPeak>, List<IPeak>> noiseFilter =
            new BinnedHighestNPeakFunction(DEFAULT_PEAKS_TO_KEEP, (int) NOISE_FILTER_INCREMENT, 0);

    /**
     * The m/z threshold to consider two peaks identical
     */
    protected final float fragmentTolerance;

    private static final PeakMzComparator peakMzComparator = new PeakMzComparator();
    private final String id;
    protected int nSpectra;
    protected float averagePrecursorMz;
    protected double sumPrecursorMz;
    protected float averagePrecursorIntens;
    protected double sumPrecursorIntens;
    protected int averageCharge;
    protected int sumCharge;
    protected ISpectrum consensusSpectrum;
    protected final List<SpectrumHolderListener> listeners = new ArrayList<SpectrumHolderListener>();

    private boolean isDirty = true;

    protected final String methodName = "Greedy Consensus Spectrum Builder";
    protected final String methodVersion = "0.1";

    /**
     * Peaks of the actual consensusSpectrum
     */
    private final List<IPeak> consensusPeaks = new ArrayList<IPeak>();

    public static final ConsensusSpectrumFactory FACTORY = new ConsensusSpectrumFactory();

    public static ConcensusSpectrumBuilderFactory buildFactory() {
        return new ConsensusSpectrumFactory();
    }

    /**
     * always use the factory to get an instance
     */
    public static class ConsensusSpectrumFactory implements ConcensusSpectrumBuilderFactory {
        private ConsensusSpectrumFactory() {
        }

        /**
         * build a new instance of the spectrum builder
         *
         * @return !null instance
         */
        @Override
        public IConsensusSpectrumBuilder getConsensusSpectrumBuilder() {
            return new GreedyConsensusSpectrum(Defaults.getFragmentIonTolerance());
        }

        public GreedyConsensusSpectrum getGreedyConsensusSpectrumBuilder() {
            return new GreedyConsensusSpectrum(Defaults.getFragmentIonTolerance());
        }

        public IConsensusSpectrumBuilder getConsensusSpectrumBuilder(String id) {
            return new GreedyConsensusSpectrum(Defaults.getFragmentIonTolerance(), id);
        }

        public GreedyConsensusSpectrum getGreedyConsensusSpectrumBuilder(String id) {
            return new GreedyConsensusSpectrum(Defaults.getFragmentIonTolerance(), id);
        }
    }

    /**
     * private to force use of the factory
     */
    private GreedyConsensusSpectrum(float fragmentTolerance) {
        this(fragmentTolerance, null);
    }

    /**
     * private to force use of the factory
     */
    private GreedyConsensusSpectrum(float fragmentTolerance, String id) {
        this.fragmentTolerance = fragmentTolerance;
        this.id = id;
    }

    public GreedyConsensusSpectrum(float fragmentTolerance, String id, int nSpectra, double sumPrecursorMz, double sumPrecursorIntens, int sumCharge, List<IPeak> peaks) {
        this.fragmentTolerance = fragmentTolerance;
        this.id = id;
        this.nSpectra = nSpectra;
        this.sumPrecursorMz = sumPrecursorMz;
        this.sumPrecursorIntens = sumPrecursorIntens;
        this.sumCharge = sumCharge;

        // update properties charge, precursor m/z and precursor intensity
        updateProperties();

        this.consensusPeaks.addAll(peaks);

        setIsDirty(true);
    }

    @Override  // TODO JG this class only correctly supports normalized spectra. Make sure the spectra are normalized
    public void addSpectra(ISpectrum... newSpectra) {
        if (newSpectra.length < 1)
            return;

        for (ISpectrum spectrum : newSpectra) {
            List<IPeak> spectrumPeaks = spectrum.getPeaks();
            addPeaksToConsensus(spectrumPeaks); // peaks are added but not additional transformation is done

            // merge identical peaks
            List<IPeak> mergedPeaks = mergeIdenticalPeaks(consensusPeaks);
            consensusPeaks.clear();
            consensusPeaks.addAll(mergedPeaks);

            sumCharge += spectrum.getPrecursorCharge();
            sumPrecursorMz += spectrum.getPrecursorMz();
            sumPrecursorIntens += 0;

            nSpectra++;
        }

        // update properties charge, precursor m/z and precursor intensity
        updateProperties();

        setIsDirty(true);

        for (SpectrumHolderListener listener : listeners)
            listener.onSpectraAdd(this, newSpectra);
    }

    public void addConsensusSpectrum(IConsensusSpectrumBuilder consensusSpectrumToAdd) {
        if (consensusSpectrumToAdd == null || consensusSpectrumToAdd.getSpectraCount() < 1)
            return;

        // add the peaks like in a "normal" spectrum - the peak count's are preserved
        addPeaksToConsensus(consensusSpectrumToAdd.getConsensusSpectrum().getPeaks());

        // merge identical peaks
        List<IPeak> mergedPeaks = mergeIdenticalPeaks(consensusPeaks);
        consensusPeaks.clear();
        consensusPeaks.addAll(mergedPeaks);

        // update the general properties
        sumCharge += consensusSpectrumToAdd.getSumCharge();
        sumPrecursorMz += consensusSpectrumToAdd.getSumPrecursorMz();
        sumPrecursorIntens += consensusSpectrumToAdd.getSumPrecursorIntensity();
        nSpectra += consensusSpectrumToAdd.getSpectraCount();

        // update properties charge, precursor m/z and precursor intensity
        updateProperties();

        setIsDirty(true);

        // this is not working correctly
        for (SpectrumHolderListener listener : listeners)
            listener.onSpectraAdd(this, consensusSpectrumToAdd.getConsensusSpectrum());
    }

    protected void updateConsensusSpectrum() {
        if (isDirty()) {

            // update the actual consensus spectrum
            List<IPeak> processedConsensusPeaks = findConsensusPeaks(consensusPeaks, nSpectra);
            consensusSpectrum = new Spectrum(id, averageCharge, averagePrecursorMz, Defaults.getDefaultQualityScorer(), processedConsensusPeaks);
            setIsDirty(false);
        }
    }

    @Override
    public void removeSpectra(ISpectrum... removed) {
        throw new UnsupportedOperationException("GreedyConsensusSpectrum does not support removing of spectra.");
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
    public void addSpectrumHolderListener(SpectrumHolderListener added) {
        listeners.add(added);
    }

    @Override
    public void removeSpectrumHolderListener(SpectrumHolderListener removed) {
        // not supported
    }

    /**
     * Adds the passed peaks to the spectrum.
     *
     * @param peaksToAdd
     */
    protected void addPeaksToConsensus(List<IPeak> peaksToAdd) {
        int posAllPeaks = 0;
        List<IPeak> newPeaks = new ArrayList<IPeak>(); // peaks with m/z values that do not yet exist

        for (int i = 0; i < peaksToAdd.size(); i++) {
            IPeak peakToAdd = peaksToAdd.get(i);
            float mzToAdd = peakToAdd.getMz();
            boolean wasAdded = false;

            for (int j = posAllPeaks; j < consensusPeaks.size(); j++) {
                IPeak currentExistingPeak = consensusPeaks.get(j);

                if (mzToAdd < currentExistingPeak.getMz()) {
                    newPeaks.add(new Peak(mzToAdd, peakToAdd.getIntensity(), peakToAdd.getCount()));
                    posAllPeaks = j;
                    wasAdded = true;
                    break;
                }

                if (mzToAdd == currentExistingPeak.getMz()) {
                    consensusPeaks.set(j, new Peak(
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
                newPeaks.add(new Peak(mzToAdd, peakToAdd.getIntensity(), peakToAdd.getCount()));
        }

        // add all new peaks
        consensusPeaks.addAll(newPeaks);
        Collections.sort(consensusPeaks, new PeakMzComparator());
    }

    /**
     * Updates all properties of the consensus spectrum as well as the actual consensus
     * spectrum.
     */
    protected void updateProperties() {
        if (nSpectra > 0) {
            averagePrecursorMz = (float) sumPrecursorMz / nSpectra;
            averageCharge = sumCharge / nSpectra;
            averagePrecursorIntens = (float) sumPrecursorIntens / nSpectra;
        } else {
            averagePrecursorMz = 0;
            averageCharge = 0;
            averagePrecursorIntens = 0;
        }
    }

    /**
     * Calls the required functions to do normalization and
     * noise filtering that result in the actual consensus
     * spectrum.
     *
     * @param input !null set of all peaks
     * @return !null set of  consensus peaks
     */
    protected static List<IPeak> findConsensusPeaks(List<IPeak> input, int nSpectra) {
        if (input.size() < 1)
            return input;

        // Step 2: adapt the peak intensities based on the probability that the peak has been observed
        List<IPeak> ret = adaptPeakIntensities(input, nSpectra);

        // Step 3: filter the spectrum
        ret = filterNoise(ret);

        return ret;
    }

    /**
     * Filters the consensus spectrum keeping only the top 5 peaks per 100 m/z
     */
    protected static List<IPeak> filterNoise(List<IPeak> inp) {
        if (inp.size() < Defaults.getDefaultConsensusMinPeaks()) {
            return inp;
        }

        // under certain conditions (averaging m/z values) the order of peaks can be disrupted
        Collections.sort(inp, peakMzComparator);
        List<IPeak> filteredSpectrum = noiseFilter.apply(inp);

        return filteredSpectrum;
    }

    /**
     * Adapt the peak intensities in consensusPeaks using the following formula:
     * I = I * (0.95 + 0.05 * (1 + pi))^5
     * where pi is the peaks probability
     */
    protected static List<IPeak> adaptPeakIntensities(List<IPeak> inp, int nSpectra) {
        List<IPeak> ret = new ArrayList<IPeak>(inp);

        for (int i = 0; i < ret.size(); i++) {
            IPeak peak = ret.get(i);
            float peakProbability = (float) peak.getCount() / (float) nSpectra;
            float newIntensity = (float) (peak.getIntensity() * (0.95 + 0.05 * Math.pow(1 + peakProbability, 5)));

            ret.set(i, new Peak(peak.getMz(), newIntensity, peak.getCount()));
        }

        return ret;
    }

    /**
     * Merges identical peaks in the consensusPeaks List based on FINAL_MZ_THRESHOLD and
     * MZ_THRESHOLD_STEP.
     */
    protected List<IPeak> mergeIdenticalPeaks(List<IPeak> inPeaks) {
        List<IPeak> filteredPeaks = new ArrayList<IPeak>();
        if (inPeaks.size() == 0)
            return filteredPeaks; // should never happen

        filteredPeaks.addAll(inPeaks);
        float mzThresholdStep = fragmentTolerance / 5; // use 4 rounds to reach the final mz threshold

        for (float range = mzThresholdStep; range < fragmentTolerance; range += mzThresholdStep) {
            List<IPeak> newPeakList = new ArrayList<IPeak>();
            IPeak currentPeak = filteredPeaks.get(0);

            for (int i = 1; i < filteredPeaks.size(); i++) {
                IPeak nextPeak = filteredPeaks.get(i);

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

            filteredPeaks.clear();
            filteredPeaks.addAll(newPeakList);
        }

        return filteredPeaks;
    }

    @Override
    public ISpectrum getConsensusSpectrum() {
        updateConsensusSpectrum();
        return consensusSpectrum;
    }


    @Override
    public void clear() {
        sumCharge = 0;
        sumPrecursorMz = 0;
        sumPrecursorIntens = 0;
        nSpectra = 0;

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

    protected boolean isDirty() {
        return isDirty;
    }

    protected void setIsDirty(boolean isDirty) {
        this.isDirty = isDirty;
    }

    @Override
    public int getSumCharge() {
        return sumCharge;
    }

    @Override
    public double getSumPrecursorMz() {
        return sumPrecursorMz;
    }

    @Override
    public double getSumPrecursorIntensity() {
        return sumPrecursorIntens;
    }

    @Override
    public List<IPeak> getRawConsensusPeaks() {
        return Collections.unmodifiableList(consensusPeaks);
    }

    @Override
    public float getFragmentIonTolerance() {
        return fragmentTolerance;
    }
}

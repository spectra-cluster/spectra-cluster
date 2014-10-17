package uk.ac.ebi.pride.spectracluster.engine;

import uk.ac.ebi.pride.spectracluster.cluster.ICluster;
import uk.ac.ebi.pride.spectracluster.similarity.ISimilarityChecker;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.util.ClusterUtilities;
import uk.ac.ebi.pride.spectracluster.util.Defaults;
import uk.ac.ebi.pride.spectracluster.util.comparator.QualityClusterComparator;

import java.util.*;

/**
 * uk.ac.ebi.pride.spectracluster.engine.PeakMatchClusteringEngine
 * performs clustering by looking at major peaks then merging clusters
 * - the this version tracks spectra already clustered
 * and attempts to combine clusters
 * <p/>
 * User: Steve
 * Date: 6/28/13
 * <p/>
 * todo: development
 */
@Deprecated
public class PeakMatchClusteringEngine implements IClusteringEngine {
    // Frank et al does 5 we do 1 more

    private final ISimilarityChecker similarityChecker;
    private final Comparator<ICluster> spectrumComparator;
    private final List<ICluster> singleSpectrumClusters = new ArrayList<ICluster>();
    private final List<ICluster> currentClusters = new ArrayList<ICluster>();
    private final Set<ISpectrum> alreadyClustered = new HashSet<ISpectrum>();
    private final IClusteringEngine clusteringEngine;

    public PeakMatchClusteringEngine(final ISimilarityChecker similarityChecker,
                                     final Comparator<ICluster> spectrumComparator,
                                     final IClusteringEngine clusteringEngine) {
        this.similarityChecker = similarityChecker;
        this.spectrumComparator = spectrumComparator;
        this.clusteringEngine = clusteringEngine;
    }


    @Override
    public ISimilarityChecker getSimilarityChecker() {
        return similarityChecker;
    }

    /**
     * Get similarity threshold used
     *
     * @return
     */
    @Override
    public double getSimilarityThreshold() {
        return clusteringEngine.getSimilarityThreshold();
    }

    /**
     * add some clusters
     */
    @Override
    public void addClusters(final ICluster... pClusters) {
        //noinspection ForLoopReplaceableByForEach
        for (int i = 0; i < pClusters.length; i++) {
            ICluster cluster = pClusters[i];
            if (cluster.getClusteredSpectraCount() == 1)
                singleSpectrumClusters.addAll(Arrays.asList(cluster));
            else
                currentClusters.add(cluster);

        }
    }


    /**
     * clusters are merged in the internal collection
     *
     * @return true is  anything happened
     */
    @Override
    public boolean processClusters() {
        if (currentClusters.isEmpty()) {
            return clusterUsingPeaks();
        } else {
            return mergeAndCombineClusters();
        }
    }

    /**
     * all subsequent passes - start  with a list if clusters and single spectrum clusters
     * merge the clusters considering only those close in mz
     *
     * @return true if anything done
     */
    protected boolean mergeAndCombineClusters() {
        int startingClusterCount = currentClusters.size() + singleSpectrumClusters.size();
        List<ICluster> singleSpectra = ClusterUtilities.removeSingleSpectrumClusters(currentClusters);
        List<ICluster> mergedClusters = mergeClusters(currentClusters);
        singleSpectra = mergeClustersWithSingles(mergedClusters, singleSpectra);

        singleSpectrumClusters.clear();
        singleSpectrumClusters.addAll(singleSpectra);

        currentClusters.clear();
        currentClusters.addAll(mergedClusters);


        int endClusterCount = currentClusters.size() + singleSpectrumClusters.size();
        return endClusterCount != startingClusterCount;
    }

    /**
     * clusters differing my more than this mx value will not be merged
     */
    public static final double MAXIMUM_MERGE_MZ_DIFFERENCE = 0.3;

    /**
     * spectra differing from clusters differing my more than this mx value will not be merged
     */
    public static final double MAXIMUM_SINGLE_SPECTRUM_MERGE_MZ_DIFFERENCE = 0.6;

    /**
     * @param mergedClusters
     * @param singleSpectra
     * @return
     */
    protected List<ICluster> mergeClustersWithSingles(List<ICluster> mergedClusters, List<ICluster> singleSpectra) {
        // let a shared function do all the dirty work so other engines can share code
        //noinspection UnnecessaryLocalVariable
        List<ICluster> retained = mergeClustersWithSingleSpectra(mergedClusters,
                singleSpectra, internalGetSimilarityChecker(), MAXIMUM_SINGLE_SPECTRUM_MERGE_MZ_DIFFERENCE);

        return retained;
    }

    /**
     * take a collection of clusters - presumably with clustered and another group - maybe
     * single spectrun clusters - that fact is not important but the groups sould be distinct
     * merge where possible and return the list of not merged spectra
     *
     * @param mergable          !null list of clusters - the clusters will change but not change in number
     * @param singles           !null list of single spectra
     * @param similarityChecker !null similarity
     * @param maxMZDIfference   ! maxMZ difference to look at for merging
     * @return !null list of non-merged clusters from singles
     */
    public List<ICluster> mergeClustersWithSingleSpectra(List<ICluster> mergable, List<ICluster> singles, ISimilarityChecker similarityChecker, double maxMZDIfference) {
        List<ICluster> retainedSingleSpectra = new ArrayList<ICluster>();
        int startIndex = 0;
        // clusters need to be compared with all clusters below them
        for (ICluster cluster : singles) {
            double currentMZ = cluster.getPrecursorMz();
            double minimimMergableMZ = currentMZ - maxMZDIfference;
            ICluster mergeWith = null;
            // start at the top mz cluster
            for (int index = startIndex; index > mergable.size(); index++) {
                ICluster test = mergable.get(index);
                if (index == startIndex) {
                    if (test.getPrecursorMz() < minimimMergableMZ) {
                        startIndex++;
                        continue; // try again
                    }
                }
                double distance = similarityChecker.assessSimilarity(test.getConsensusSpectrum(), cluster.getConsensusSpectrum());
                if (distance <= clusteringEngine.getSimilarityThreshold()) {
                    mergeWith = test;
                    break; // found who to merge with
                }

            }
            if (mergeWith == null) {
                retainedSingleSpectra.add(cluster); // nothing to merge with so keep single
            } else {  // merge with a close enough cluster
                // note this may disturb the order of the cluster list but should not stop
                // the algorithm from working
                final List<ISpectrum> clusteredSpectra = cluster.getClusteredSpectra();
                mergeWith.addSpectra(clusteredSpectra.toArray(new ISpectrum[clusteredSpectra.size()]));
            }
        }


        // make sure the mergable are still in mz order
        Collections.sort(mergable);


        // make sure the retainedSingleSpectra are still in mz order
        Collections.sort(retainedSingleSpectra);
        return retainedSingleSpectra;
    }


    /**
     * merge clusters among themselves returning a list of surviving clusters which will replace mergable
     *
     * @param mergable
     * @return
     */
    protected List<ICluster> mergeClusters(List<ICluster> mergable) {
        List<ICluster> retained = new ArrayList<ICluster>();
        // clusters need to be compared with all clusters below them
        for (ICluster cluster : mergable) {
            double currentMZ = cluster.getPrecursorMz();
            double minimimMergableMZ = currentMZ - MAXIMUM_MERGE_MZ_DIFFERENCE;
            ICluster mergeWith = null;
            if (retained.isEmpty()) {   // start with the first cluster
                retained.add(cluster);
                continue;
            }

            // start at the top mz cluster
            for (int index = retained.size() - 1; index >= 0; index--) {
                ICluster test = retained.get(index);
                if (test.getPrecursorMz() < minimimMergableMZ)
                    break; // no more to consider
                final ISpectrum cs1 = test.getConsensusSpectrum();
                double similarity = similarityChecker.assessSimilarity(cs1, cluster.getConsensusSpectrum());
                if (similarity >= Defaults.getSimilarityThreshold()) {
                    mergeWith = test;
                    break; // found who to merge with
                }

            }
            if (mergeWith == null) {
                retained.add(cluster); // nothing to merge with so keep the cluster
            } else {  // merge with a close enough cluster
                final List<ISpectrum> clusteredSpectra = cluster.getClusteredSpectra();
                mergeWith.addSpectra(clusteredSpectra.toArray(new ISpectrum[clusteredSpectra.size()]));
            }
        }
        // make sure the retained are still in mz order
        Collections.sort(retained);
        return retained;
    }

    /**
     * generate clusters in current Clusters from read clusters
     *
     * @return always true something is dome
     */
    protected boolean clusterUsingPeaks() {
        Collections.sort(singleSpectrumClusters, QualityClusterComparator.INSTANCE);   // sort by quality
        for (int index = 0; index < singleSpectrumClusters.size(); index++) {
            ICluster readCluster = singleSpectrumClusters.get(index);
            if (readCluster.getClusteredSpectraCount() != 1)
                throw new IllegalStateException("this should be a a single spectrum cluster"); // ToDo change
            final ISpectrum theSpectrum = readCluster.getHighestQualitySpectrum();
            final int[] peaks = theSpectrum.asMajorPeakMZs(Defaults.getMajorPeakCount());
            //noinspection ForLoopReplaceableByForEach
            for (int i = 0; i < peaks.length; i++) {
                if (alreadyClustered.contains(theSpectrum))    // we are already in a cluster
                    continue;
                int peak = peaks[i];
                for (int index2 = index; index2 < singleSpectrumClusters.size(); index2++) {
                    ICluster addedCluster = singleSpectrumClusters.get(index2);
                    ISpectrum addedSpectrum = readCluster.getHighestQualitySpectrum();
                    if (alreadyClustered.contains(addedSpectrum))
                        continue;
                    if (!addedSpectrum.containsMajorPeak(peak, Defaults.getMajorPeakCount()))   // we do not have the peak
                        continue;
                    clusteringEngine.addClusters(addedCluster);
                }
                if (clusteringEngine.size() < 2)
                    continue; // nothing to cluster
                clusteringEngine.processClusters();
                final Collection<ICluster> clusters = clusteringEngine.getClusters();
                currentClusters.addAll(clusters);
                for (ICluster cluster : clusters) {
                    alreadyClustered.addAll(cluster.getClusteredSpectra());
                }
            }

        }

        alreadyClustered.clear(); // we don't need to remember these
        return true;
    }


    /**
     * Get clustered clusters
     * SLewis - I think a guarantee that they are sorted by MZ is useful
     */
    @Override
    public List<ICluster> getClusters() {
        if (currentClusters.isEmpty()) {        // pass 1
            List<ICluster> ret = new ArrayList<ICluster>(singleSpectrumClusters);
            Collections.sort(ret);
            return ret;
        } else {
            List<ICluster> ret = new ArrayList<ICluster>(currentClusters);
            ret.addAll(singleSpectrumClusters);
            Collections.sort(ret);
            return ret;

        }
    }


    /**
     * used to expose internals for overriding classes only
     *
     * @return
     */
    protected ISimilarityChecker internalGetSimilarityChecker() {
        return similarityChecker;
    }

    /**
     * used to expose internals for overriding classes only
     *
     * @return
     */
    @SuppressWarnings("UnusedDeclaration")
    protected Comparator<ICluster> internalGetSpectrumComparator() {
        return spectrumComparator;
    }

    /**
     * allow engines to be named
     *
     * @return
     */
    @Override
    public String toString() {
        String name = this.getClass().getName();
        if (name != null)
            return name;
        return super.toString();
    }

    /**
     * total number of clusters including queued clustersToAdd
     *
     * @return
     */
    @Override
    public int size() {
        return singleSpectrumClusters.size();  // todo do better
    }
}

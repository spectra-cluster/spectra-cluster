package uk.ac.ebi.pride.spectracluster.engine;

import uk.ac.ebi.pride.spectracluster.cluster.ICluster;
import uk.ac.ebi.pride.spectracluster.similarity.ISimilarityChecker;

import java.util.*;

/**
 * uk.ac.ebi.pride.spectracluster.engine.WrappedIncrementalClusteringEngine
 * Wraps an  IIncrementalClusteringEngine to do one pass clustering
 * User: Steve
 * Date: 8/13/13
 * <p/>
 * todo: development, since it is never used
 */
@Deprecated
public class WrappedIncrementalClusteringEngine implements IClusteringEngine {

    private boolean dirty;
    private final IIncrementalClusteringEngine realEngine;
    private final List<ICluster> clusters = new ArrayList<ICluster>();

    public WrappedIncrementalClusteringEngine(final IIncrementalClusteringEngine pRealEngine) {
        realEngine = pRealEngine;
    }

    public IIncrementalClusteringEngine getRealEngine() {
        return realEngine;
    }

    public boolean isDirty() {
        return dirty;
    }

    public void setDirty(final boolean pDirty) {
        dirty = pDirty;
    }

    /**
     * simple get of raw clusters array for internal use
     *
     * @return
     */
    protected List<ICluster> internalGetClusters() {
        return clusters;
    }

    /**
     * Get clustered clusters sorted by MZ is useful
     *
     * @return !null list this will be sorted by mz a include clusters of all sizes
     */
    @Override
    public List<ICluster> getClusters() {
        final List<ICluster> internalClusters = internalGetClusters();
        Set<ICluster> internalSet = new HashSet<ICluster>(internalClusters);
        IIncrementalClusteringEngine engine = getRealEngine();
        final Collection<ICluster> lastClusters = engine.getClusters();
        internalSet.addAll(lastClusters);
        List<ICluster> ret = new ArrayList<ICluster>(internalSet);
        Collections.sort(ret);
        return ret;
    }

    /**
     * add some clusters
     */
    @Override
    public void addClusters(final ICluster... cluster) {
        final List<ICluster> internalClusters = internalGetClusters();
        IIncrementalClusteringEngine engine = getRealEngine();
        //noinspection ForLoopReplaceableByForEach
        for (int i = 0; i < cluster.length; i++) {
            ICluster added = cluster[i];
            final Collection<ICluster> finalClusters = engine.addClusterIncremental(added);
            if (!finalClusters.isEmpty())
                internalClusters.addAll(finalClusters);
        }
        setDirty(true);
    }

    /**
     * clusters are merged in the internal collection
     * processing is incremental
     *
     * @return true is  anything happened
     */
    @Override
    public boolean processClusters() {
        boolean wasDirty = isDirty();
        setDirty(false);
        return wasDirty;
    }

    @Override
    public ISimilarityChecker getSimilarityChecker() {
        return realEngine.getSimilarityChecker();
    }

    /**
     * Get similarity threshold used
     *
     * @return
     */
    @Override
    public double getSimilarityThreshold() {
        return realEngine.getSimilarityThreshold();
    }

    /**
     * total number of clusters including queued clustersToAdd
     *
     * @return
     */
    @Override
    public int size() {
        return internalGetClusters().size();
    }
}

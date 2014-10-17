package uk.ac.ebi.pride.spectracluster.engine;

/**
 * Factory for making WrappedIncrementalClusteringEngine
 * <p/>
 * The only reason to keep this factory class is for the default configuration
 *
 * @author Steve Lewis
 * @author Rui Wang
 * @version $Id$
 *          <p/>
 *          todo: development since it is only used in unit tests
 */
@Deprecated
public class WrappedIncrementalClusteringEngineFactory {

    /**
     * make a copy of the clustering engine
     *
     * @return
     */
    public IClusteringEngine getClusteringEngine(Object... otherdata) {
        if (otherdata.length < 1)
            throw new IllegalArgumentException("WrappedClusteringEngine needs a Double as WindowSize"); //
        double ws = (Double) otherdata[0];
        float windowSize = (float) ws;
        final IIncrementalClusteringEngine incrementalClusteringEngine = EngineFactories.buildIncrementalClusteringEngineFactory(windowSize).buildInstance();
        return new WrappedIncrementalClusteringEngine(incrementalClusteringEngine);
    }
}

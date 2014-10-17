package uk.ac.ebi.pride.spectracluster.engine;

import uk.ac.ebi.pride.spectracluster.similarity.ISimilarityChecker;
import uk.ac.ebi.pride.spectracluster.util.Defaults;
import uk.ac.ebi.pride.spectracluster.util.IDefaultingFactory;
import uk.ac.ebi.pride.spectracluster.util.comparator.ClusterComparator;

/**
 * Factory for making PeakMatchClusteringEngine
 * <p/>
 * The only reason to keep this factory class is for the default configuration
 *
 * @author Rui Wang
 * @version $Id$
 *          <p/>
 *          todo: development
 */
@Deprecated
public class PeakMatchClusteringEngineFactory {

    /**
     * build a new version
     *
     * @return
     */
    public PeakMatchClusteringEngine getPeakMatchClusteringEngine() {
        final ISimilarityChecker similarityChecker = Defaults.getDefaultSimilarityChecker();
        final ClusterComparator comparator = Defaults.getDefaultSpectrumComparator();
        final IDefaultingFactory<IClusteringEngine> engineIDefaultingFactory = EngineFactories.buildClusteringEngineFactory(similarityChecker, comparator, Defaults.getSimilarityThreshold());
        final IClusteringEngine clusteringEngine = engineIDefaultingFactory.buildInstance();
        return new PeakMatchClusteringEngine(similarityChecker, comparator, clusteringEngine);
    }
}

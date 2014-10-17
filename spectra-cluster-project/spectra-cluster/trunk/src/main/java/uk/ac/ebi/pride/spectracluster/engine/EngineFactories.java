package uk.ac.ebi.pride.spectracluster.engine;

import uk.ac.ebi.pride.spectracluster.similarity.ISimilarityChecker;
import uk.ac.ebi.pride.spectracluster.util.Defaults;
import uk.ac.ebi.pride.spectracluster.util.IDefaultingFactory;
import uk.ac.ebi.pride.spectracluster.util.comparator.ClusterComparator;

/**
 * uk.ac.ebi.pride.spectracluster.engine.EngineFactories
 *
 * @author Steve Lewis
 * @date 28/05/2014
 */
public class EngineFactories {

    /**
     * make a clustering engine assuming the defaults are used
     */
    public static IDefaultingFactory<IClusteringEngine> DEFAULT_CLUSTERING_ENGINE_FACTORY =
            new ClusteringEngineFactory(Defaults.getDefaultSimilarityChecker(),
                    Defaults.getDefaultSpectrumComparator(),
                    Defaults.getSimilarityThreshold());


    public static IDefaultingFactory<IClusteringEngine> buildClusteringEngineFactory(
            final ISimilarityChecker pSimilarityChecker,
            ClusterComparator pSpectrumComparator,
            double threshold) {
        return new ClusteringEngineFactory(pSimilarityChecker,
                pSpectrumComparator,
                threshold);
    }


    protected static class ClusteringEngineFactory implements IDefaultingFactory<IClusteringEngine> {
        private final ISimilarityChecker similarityChecker;
        private final ClusterComparator spectrumComparator;
        private final double similarityThreshold;

        private ClusteringEngineFactory(final ISimilarityChecker pSimilarityChecker,
                                        ClusterComparator pSpectrumComparator,
                                        double threshold
        ) {
            similarityChecker = pSimilarityChecker;
            spectrumComparator = pSpectrumComparator;
            similarityThreshold = threshold;
        }

        /**
         * make a copy of the clustering engine
         *
         * @return
         */
        @Override
        public IClusteringEngine buildInstance(Object... otherdata) {
            return new ClusteringEngine(similarityChecker, spectrumComparator, similarityThreshold);
        }
    }

    public static IDefaultingFactory<IClusteringEngine> DEFAULT_INCREMENTAL_CLUSTERING_ENGINE_FACTORY =
            new ClusteringEngineFactory(Defaults.getDefaultSimilarityChecker(),
                    Defaults.getDefaultSpectrumComparator(),
                    Defaults.getSimilarityThreshold());

    /**
     * make a clustering engine assuming the defaults are used  for all but window size
     */
    public static IDefaultingFactory<IIncrementalClusteringEngine> buildIncrementalClusteringEngineFactory(
            float windowSize) {
        return new IncrementalClusteringEngineFactory(Defaults.getDefaultSimilarityChecker(),
                Defaults.getDefaultSpectrumComparator(),
                Defaults.getSimilarityThreshold(),
                windowSize);
    }


    public static IDefaultingFactory<IIncrementalClusteringEngine> buildIncrementalClusteringEngineFactory(
            final ISimilarityChecker pSimilarityChecker,
            ClusterComparator pSpectrumComparator,
            double threshold,
            float windowSize) {
        return new IncrementalClusteringEngineFactory(pSimilarityChecker,
                pSpectrumComparator,
                threshold,
                windowSize);
    }


    public static class IncrementalClusteringEngineFactory implements IDefaultingFactory<IIncrementalClusteringEngine> {
        private final ISimilarityChecker similarityChecker;
        private final ClusterComparator spectrumComparator;
        private final double similarityThreshold;
        private final float windowSize;

        private IncrementalClusteringEngineFactory(final ISimilarityChecker pSimilarityChecker,
                                                   ClusterComparator pSpectrumComparator,
                                                   double threshold,
                                                   float windowSize
        ) {
            similarityChecker = pSimilarityChecker;
            spectrumComparator = pSpectrumComparator;
            similarityThreshold = threshold;
            this.windowSize = windowSize;
        }

        /**
         * make a copy of the clustering engine
         *
         * @return
         */
        public IIncrementalClusteringEngine getIncrementalClusteringEngine(float ws) {
            return new IncrementalClusteringEngine(similarityChecker, spectrumComparator, ws, similarityThreshold);
        }

        /**
         * make a copy of the clustering engine
         *
         * @return
         */
        @Override
        public IIncrementalClusteringEngine buildInstance(Object... otherdata) {
            return new IncrementalClusteringEngine(similarityChecker, spectrumComparator, windowSize, similarityThreshold);
        }
    }

}

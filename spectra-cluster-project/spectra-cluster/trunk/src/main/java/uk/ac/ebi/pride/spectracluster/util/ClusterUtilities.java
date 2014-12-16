package uk.ac.ebi.pride.spectracluster.util;

import uk.ac.ebi.pride.spectracluster.cluster.ICluster;
import uk.ac.ebi.pride.spectracluster.cluster.SpectralCluster;
import uk.ac.ebi.pride.spectracluster.similarity.ISimilarityChecker;
import uk.ac.ebi.pride.spectracluster.spectrum.*;

import javax.annotation.Nonnull;
import javax.annotation.Nullable;
import java.util.*;


/**
 * uk.ac.ebi.pride.spectracluster.util.ClusterUtilities
 * a list of stateless static functions for manipulating clusters, lists of clusters
 * and performing other common chores
 *
 * @author Steve Lewis
 * @date 5/10/13
 */
public final class ClusterUtilities {

    private ClusterUtilities() {
    }


    /**
     * return non-null if c1 contains all spactra in c2 or vice
     * versa
     *
     * @param c1 cluster
     * @param c2 cluster
     * @return non-null if we can use the return as an enclosing cluster
     */
    public static
    @Nullable
    ICluster clusterFullyContains(@Nonnull ICluster c1, @Nonnull ICluster c2) {
        int size1 = c1.getClusteredSpectraCount();
        int size2 = c2.getClusteredSpectraCount();
        if (size1 == size2) {
            if (c1.getSpectralId().equals(c2.getSpectralId()))
                return c1; // same size same spectra;
        }
        ICluster containing = c1;
        ICluster smaller = c2;
        if (size1 < size2) {
            containing = c2;
            smaller = c1;
        }
        Set<String> spectralIds = containing.getSpectralIds();
        if (spectralIds.containsAll(smaller.getSpectralIds()))
            return containing;   // contains all spectra in c1

        return null; // clusters are distinct
    }


    /**
     * score commonality
     *
     * @param existing cluster
     * @param added    cluster
     * @return non-null if we can use the return as an enclosing cluster
     */
    @Deprecated // TODO JG function highly similar to IncrementalClusteringEngine::getSharedSpectraIds
    public static double clusterFullyContainsScore(@Nonnull ICluster existing, @Nonnull ICluster added) {
        Set<String> spectralIds1 = existing.getSpectralIds();
        Set<String> spectralIds2 = added.getSpectralIds();
        double minSize = Math.min(spectralIds1.size(), spectralIds2.size());

        Set<String> common = new HashSet<String>(spectralIds1);
        common.retainAll(spectralIds2);

        return (double) common.size() / minSize;
    }

    /**
     * expose critical code for demerge - THIS NEVER CHANGES INTERNAL STATE and
     * usually is called on removed clusters
     *
     * @return !null Cluster
     */
    public static List<ICluster> findNoneFittingSpectra(ICluster cluster, ISimilarityChecker similarityChecker, double threshold) {
        List<ICluster> noneFittingSpectra = new ArrayList<ICluster>();

        if (cluster.getClusteredSpectra().size() > 1) {
            for (ISpectrum spectrum : cluster.getClusteredSpectra()) {
                final ISpectrum consensusSpectrum = cluster.getConsensusSpectrum();
                final double similarityScore = similarityChecker.assessSimilarity(consensusSpectrum, spectrum);

                if (similarityScore < threshold) {
                    noneFittingSpectra.add(ClusterUtilities.asCluster(spectrum));
                }
            }
        }

        return noneFittingSpectra;
    }

    /**
     * allow nonfitting spectra to leave and retuen a list of clusters to write out
     *
     * @param cluster
     * @return !null List<ISpectralCluster
     */
    @Nonnull
    public static List<ICluster> removeNonFittingSpectra(@Nonnull ICluster cluster, @Nonnull ISimilarityChecker similarityChecker, double threshold) {
        final List<ICluster> allClusters = findNoneFittingSpectra(cluster, similarityChecker, threshold);
        List<ICluster> holder = new ArrayList<ICluster>();
        if (!allClusters.isEmpty()) {
            for (ICluster removedCluster : allClusters) {

                // drop all spectra
                final List<ISpectrum> clusteredSpectra = removedCluster.getClusteredSpectra();
                ISpectrum[] allRemoved = clusteredSpectra.toArray(new ISpectrum[clusteredSpectra.size()]);
                cluster.removeSpectra(allRemoved);

                // and save as stand alone
                holder.add(removedCluster);
            }

        }
        if (cluster.getClusteredSpectraCount() > 0)
            holder.add(cluster);
        return holder;
    }


    /**
     * find the highest quality spectrum in a list of clusters
     *
     * @param copied should be non-empty array
     * @return !null spectrum unless copied is empty
     */
    public static ISpectrum getHighestQualitySpectrum(ICluster... copied) {
        if (copied.length == 0)
            return null;
        ISpectrum ret = copied[0].getHighestQualitySpectrum();
        for (int i = 1; i < copied.length; i++) {
            ISpectrum highestQualitySpectrum = copied[i].getHighestQualitySpectrum();
            if (!ret.equivalent(highestQualitySpectrum))
                throw new IllegalStateException("AlternativeSpectralClusters MUST have the same highest quality spectrum");
        }

        return ret;
    }


    /**
     * take out all clusters consisting of a single spectrum and return them as a list
     * - engines will do this as a step in reclustering
     *
     * @param clusters !null   a list of clusters - this WILL be modified
     * @return !null list of clusters containing a single spectrum
     */
    public static List<ICluster> removeSingleSpectrumClusters(List<ICluster> clusters) {
        return removeSingleSpectrumClustersSizedLessThan(clusters, 1);
    }

    /**
     * take out all clusters sized les shtne size and return them as a list of single spectrum clusters
     * - engines will do this as a step in reclustering
     *
     * @param pClusters !null   a list of clusters - this WILL be modified
     * @param size      - clusters sized less than or equal to this will be removed and returned as single spectrum clusters
     * @return !null list of clusters containing a single spectrum
     */
    public static List<ICluster> removeSingleSpectrumClustersSizedLessThan(final List<ICluster> pClusters, final int size) {
        List<ICluster> retained = new ArrayList<ICluster>();
        List<ICluster> asSingleSpectra = new ArrayList<ICluster>();
        for (ICluster cluster : pClusters) {
            if (cluster.getClusteredSpectraCount() <= size) {
                final List<ISpectrum> clusteredSpectra = cluster.getClusteredSpectra();
                for (ISpectrum spectrum : clusteredSpectra) {
                    asSingleSpectra.add(ClusterUtilities.asCluster(spectrum)); // di not retain but return as one spectrum clusters
                }
            } else {
                retained.add(cluster); // large enough keep it
            }
        }
        pClusters.clear();
        pClusters.addAll(retained);

        return asSingleSpectra;
    }

    /**
     * convert a spectrum into cluster
     *
     * @param spectrum given spectrum
     * @return a spectral cluster
     */
    public static ICluster asCluster(ISpectrum spectrum) {
        ICluster ret = new SpectralCluster((String) null, Defaults.getDefaultConsensusSpectrumBuilder());
        ret.addSpectra(spectrum);
        return ret;
    }


    /**
     * return a list of all spectra in the list of spectra sorted by charge then mz
     *
     * @param spectra !null list of spectra
     * @return !null list of spectra
     */
    public static List<ICluster> asClusters(List<ISpectrum> spectra) {
        List<ICluster> holder = new ArrayList<ICluster>();
        for (ISpectrum spectrum : spectra) {
            holder.add(asCluster(spectrum));
        }
        Collections.sort(holder);
        return holder;
    }


    /**
     * return a comma delimited set of most common peptides
     *
     * @return as above
     */
    public static String mostCommonPeptides(ICluster cluster) {
        String property = cluster.getProperty(KnownProperties.MOST_COMMON_PEPTIDE_KEY);
        if(property != null)
            return property;

        final List<ISpectrum> clusteredSpectra = cluster.getClusteredSpectra();
        //noinspection UnnecessaryLocalVariable
        String petides = SpectrumUtilities.mostCommonPeptides(clusteredSpectra);
        return petides;
    }

//         todo to we need this 4-jun-2014
//    public static String[] getMostCommonPeptides(ICluster cluster) {
//        final List<ISpectrum> clusteredSpectra = cluster.getClusteredSpectra();
//        //noinspection UnnecessaryLocalVariable
//        final List<String> peptideList = SpectrumUtilities.getPeptideList(clusteredSpectra);
//        //noinspection UnnecessaryLocalVariable,UnusedDeclaration,UnusedAssignment
//        final String[] stringsByOccurance = CountedString.getStringsByOccurance(peptideList);
//
//        return stringsByOccurance;
//    }


    /**
     * return a list of all spectra in the list of clusters sorted by charge then mz
     *
     * @param clusters !null list of clusters
     * @return !null list of spectra
     */
    public static List<ISpectrum> extractSpectra(List<ICluster> clusters) {
        List<ISpectrum> holder = new ArrayList<ISpectrum>();
        for (ICluster cluster : clusters) {
            holder.addAll(cluster.getClusteredSpectra());
        }
        Collections.sort(holder);
        return holder;
    }
}

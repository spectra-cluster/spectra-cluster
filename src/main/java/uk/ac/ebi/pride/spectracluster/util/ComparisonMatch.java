package uk.ac.ebi.pride.spectracluster.util;

/**
 * This class only stores a certain similarity and the
 * spectrum id it was associated with. It is intended to
 * store similarity matches within a cluster.
 *
 * Created by jg on 06.05.15.
 */
public class ComparisonMatch implements Comparable<ComparisonMatch> {
    private final String spectrumId;
    /**
     * Single precision is sufficient for this
     */
    private final float similarity;

    public ComparisonMatch(String spectrumId, float similarity) {
        this.spectrumId = spectrumId;
        this.similarity = similarity;
    }

    public String getSpectrumId() {
        return spectrumId;
    }

    public float getSimilarity() {
        return similarity;
    }

    @Override
    public int compareTo(ComparisonMatch o) {
        return Float.compare(this.similarity, o.similarity);
    }
}

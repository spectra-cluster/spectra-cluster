package uk.ac.ebi.pride.spectracluster.util;

/**
 * Equivalent
 * like Comparable or equals but weaker - used for similarity comparisons
 * it says an object can say it is equivalent to another object - not
 * necessarily of the same Class - frequently implementing the same interface
 *
 * @author Steve Lewis
 */
public interface Equivalent<T> {

    /**
     * like equals but weaker - says other is equivalent to this
     *
     * @param other poiibly null other object
     * @return true if other is "similar enough to this"
     */
    public boolean equivalent(T other);
}

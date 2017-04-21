package uk.ac.ebi.pride.spectracluster.util;

/**
 * uk.ac.ebi.pride.spectracluster.util.IDefaultingFactory
 *
 * @author Steve Lewis
 */
public interface IDefaultingFactory<T> {
    /**
     * create an object of type T using defaults
     *
     * @return
     */
    T buildInstance(Object... input);
}

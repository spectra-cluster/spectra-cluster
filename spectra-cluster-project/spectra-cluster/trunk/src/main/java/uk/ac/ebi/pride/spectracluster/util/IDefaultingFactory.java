package uk.ac.ebi.pride.spectracluster.util;

/**
 * uk.ac.ebi.pride.spectracluster.util.IDefaultingFactory
 *
 * @author Steve Lewis
 * @date 28/05/2014
 */
public interface IDefaultingFactory<T> {
    /**
     * create an object of type T using defaults
     *
     * @return
     */
    public T buildInstance(Object... input);
}

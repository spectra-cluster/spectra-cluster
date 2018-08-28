package uk.ac.ebi.pride.spectracluster.util;

import java.io.Serializable;

/**
 * uk.ac.ebi.pride.spectracluster.util.IDefaultingFactory
 *
 * @author Steve Lewis
 */
public interface IDefaultingFactory<T> extends Serializable{
    /**
     * create an object of type T using defaults
     *
     * @return
     */
    T buildInstance(Object... input);
}

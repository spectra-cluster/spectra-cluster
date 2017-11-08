package uk.ac.ebi.pride.spectracluster.util.predicate;

import java.io.Serializable;

/**
 * Compares two spectra
 *
 * Created by jg on 20.05.15.
 * @author Yasset Perez-Riverol
 */
public interface IComparisonPredicate<T> extends Serializable{
    boolean apply(T o1, T o2);
}

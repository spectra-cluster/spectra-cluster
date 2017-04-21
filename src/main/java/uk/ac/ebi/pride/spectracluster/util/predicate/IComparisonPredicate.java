package uk.ac.ebi.pride.spectracluster.util.predicate;

/**
 * Compares two spectra
 *
 * Created by jg on 20.05.15.
 */
public interface IComparisonPredicate<T> {
    boolean apply(T o1, T o2);
}

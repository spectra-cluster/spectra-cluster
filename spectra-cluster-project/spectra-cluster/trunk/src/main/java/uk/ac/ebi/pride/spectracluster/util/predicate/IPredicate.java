package uk.ac.ebi.pride.spectracluster.util.predicate;

/**
 * Determine a true or false value for a given input
 *
 * @author Rui Wang
 * @version $Id$
 */
public interface IPredicate<T> {

    boolean apply(T o);
}

package uk.ac.ebi.pride.spectracluster.util.predicate;

import java.io.Serializable;

/**
 * Determine a true or false value for a given input
 *
 * @author Rui Wang
 * @author Yasset Perez-Riverol
 * @version $Id$
 */
public interface IPredicate<T> extends Serializable{

    boolean apply(T o);
}

package uk.ac.ebi.pride.spectracluster.util.predicate;

import java.util.Arrays;

/**
 * Utility class for working with Predicates
 *
 * @author Rui Wang
 * @version $Id$
 */
public final class ComparisonPredicates {

    /**
     * Perform AND operations on a group of predicates
     *
     * @param predicates a group of input predicates
     * @param <T>        predicate input object
     * @return a new predicate represents the AND of all predicates
     */
    public static <T> IComparisonPredicate<T> and(final Iterable<? extends IComparisonPredicate<? super T>> predicates) {
        return new IComparisonPredicate<T>() {
            @Override
            public boolean apply(T obj1, T obj2) {
                for (IComparisonPredicate<? super T> predicate : predicates) {
                    if (!predicate.apply(obj1, obj2))
                        return false;
                }
                return true;
            }
        };
    }

    /**
     * Perform AND operations on a group of predicates
     *
     * @param predicates a group of input predicates
     * @param <T>        predicate input object
     * @return a new predicate represents the AND of all predicates
     */
    public static <T> IComparisonPredicate<T> and(final IComparisonPredicate<? super T>... predicates) {
        return and(Arrays.asList(predicates));
    }

    /**
     * Perform OR operations on a group of predicates
     *
     * @param predicates a group of input predicates
     * @param <T>        predicate input object
     * @return a new predicate represents the OR of all predicates
     */
    public static <T> IComparisonPredicate<T> or(final Iterable<? extends IComparisonPredicate<? super T>> predicates) {
        return new IComparisonPredicate<T>() {
            @Override
            public boolean apply(T obj1, T obj2) {
                for (IComparisonPredicate<? super T> predicate : predicates) {
                    if (predicate.apply(obj1, obj2))
                        return true;
                }
                return false;
            }
        };
    }

    /**
     * Perform OR operations on a group of predicates
     *
     * @param predicates a group of input predicates
     * @param <T>        predicate input object
     * @return a new predicate represents the OR of all predicates
     */
    public static <T> IComparisonPredicate<T> or(final IComparisonPredicate<? super T>... predicates) {
        return or(Arrays.asList(predicates));
    }

    /**
     * Return a predicate which always true
     *
     * @param <T> input object
     * @return predicate that always returns true
     */
    public static <T> IComparisonPredicate<T> alwaysTrue() {
        return new IComparisonPredicate<T>() {
            @Override
            public boolean apply(T obj1, T obj2) {
                return true;
            }
        };
    }


    /**
     * Return a predicate which always false
     *
     * @param <T> input object
     * @return predicate that always returns false
     */
    public static <T> IComparisonPredicate<T> alwaysFalse() {
        return new IComparisonPredicate<T>() {
            @Override
            public boolean apply(T obj1, T obj2) {
                return false;
            }
        };
    }

}

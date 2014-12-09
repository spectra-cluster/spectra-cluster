package uk.ac.ebi.pride.spectracluster.util.predicate;

import uk.ac.ebi.pride.spectracluster.util.function.IFunction;

import java.util.Arrays;

/**
 * Utility class for working with Predicates
 *
 * @author Rui Wang
 * @version $Id$
 */
public final class Predicates {

    /**
     * Perform AND operations on a group of predicates
     *
     * @param predicates a group of input predicates
     * @param <T>        predicate input object
     * @return a new predicate represents the AND of all predicates
     */
    public static <T> IPredicate<T> and(final Iterable<? extends IPredicate<? super T>> predicates) {
        return new IPredicate<T>() {
            @Override
            public boolean apply(T obj) {
                for (IPredicate<? super T> predicate : predicates) {
                    if (!predicate.apply(obj))
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
    public static <T> IPredicate<T> and(final IPredicate<? super T>... predicates) {
        return and(Arrays.asList(predicates));
    }

    /**
     * Perform OR operations on a group of predicates
     *
     * @param predicates a group of input predicates
     * @param <T>        predicate input object
     * @return a new predicate represents the OR of all predicates
     */
    public static <T> IPredicate<T> or(final Iterable<? extends IPredicate<? super T>> predicates) {
        return new IPredicate<T>() {
            @Override
            public boolean apply(T obj) {
                for (IPredicate<? super T> predicate : predicates) {
                    if (predicate.apply(obj))
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
    public static <T> IPredicate<T> or(final IPredicate<? super T>... predicates) {
        return or(Arrays.asList(predicates));
    }


    /**
     * Combine a predicate and a function together
     *
     * @param predicate given predicate
     * @param func      given function to transform the object
     * @param <A>       input object
     * @param <B>       function's output object
     * @return boolean
     */
    public static <A, B> IPredicate<A> condition(final IPredicate<B> predicate, final IFunction<A, ? extends B> func) {
        return new IPredicate<A>() {
            @Override
            public boolean apply(A obj) {
                B result = func.apply(obj);
                return predicate.apply(result);
            }
        };
    }

    /**
     * Return a predicate which always true
     *
     * @param <T> input object
     * @return predicate that always returns true
     */
    public static <T> IPredicate<T> alwaysTrue() {
        return new IPredicate<T>() {
            @Override
            public boolean apply(T obj) {
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
    public static <T> IPredicate<T> alwaysFalse() {
        return new IPredicate<T>() {
            @Override
            public boolean apply(T obj) {
                return false;
            }
        };
    }

}

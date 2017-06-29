package uk.ac.ebi.pride.spectracluster.util;

import uk.ac.ebi.pride.spectracluster.util.function.IFunction;
import uk.ac.ebi.pride.spectracluster.util.predicate.IPredicate;

import java.util.Collection;
import java.util.stream.Collectors;

/**
 * Utilities class for filtering and transforming a collection
 *
 * @author Rui Wang
 * @version $Id$
 */
public final class CollectionUtils {

    /**
     * Filter an input collection using a given predicate
     *
     * @param input     input collection
     * @param predicate predicate
     * @param <T>       generic object
     */
    public static <T> void filter(Collection<T> input, IPredicate<T> predicate) {
        input.removeIf(element -> !predicate.apply(element));
    }

    /**
     * Transform an input collection to a new collection using a given function
     *
     * @param input input collection
     * @param func  function to be applied
     * @param <A>   function input element type
     * @param <B>   function output element type
     * @return a transformed collection
     */
    public static <A, B> Collection<B> transform(Collection<A> input, IFunction<A, B> func) {
        Collection<B> results = input.stream().map(func::apply).collect(Collectors.toList());
        return results;
    }

}

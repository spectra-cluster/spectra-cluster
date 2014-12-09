package uk.ac.ebi.pride.spectracluster.util.function;

import uk.ac.ebi.pride.spectracluster.util.predicate.IPredicate;

/**
 * Utility class for working with Functions
 *
 * @author Rui Wang
 * @version $Id$
 */
public final class Functions {

    /**
     * Compose two functions together
     *
     * @param startFunc start function which takes input
     * @param endFunc   end function which produce output
     * @param <A>       input object
     * @param <B>       intermediate object
     * @param <C>       output object
     * @return C as the output
     */
    public static <A, B, C> IFunction<A, C> compose(final IFunction<A, ? extends B> startFunc, final IFunction<B, C> endFunc) {
        return new IFunction<A, C>() {
            @Override
            public C apply(A obj) {
                B startFuncResult = startFunc.apply(obj);
                return endFunc.apply(startFuncResult);
            }
        };
    }

    /**
     * Conditioned applying function on object using a predicate
     * If predicate is satisfied then return the function result, otherwise,
     * return null.
     *
     * @param func      given function
     * @param predicate given predicate
     * @param <A>       input object
     * @param <B>       output object
     * @return B as the output
     */
    public static <A, B> IFunction<A, B> condition(final IFunction<A, ? extends B> func, final IPredicate<B> predicate) {
        return new IFunction<A, B>() {
            @Override
            public B apply(A obj) {
                B result = func.apply(obj);
                return predicate.apply(result) ? result : null;
            }
        };
    }
}

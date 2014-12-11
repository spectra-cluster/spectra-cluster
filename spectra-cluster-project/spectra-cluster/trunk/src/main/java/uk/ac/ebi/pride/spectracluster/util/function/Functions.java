package uk.ac.ebi.pride.spectracluster.util.function;

import uk.ac.ebi.pride.spectracluster.util.predicate.IPredicate;

import java.util.Arrays;

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
    public static <A, B, C> IFunction<A, C> join(final IFunction<A, ? extends B> startFunc, final IFunction<B, C> endFunc) {
        return new IFunction<A, C>() {
            @Override
            public C apply(A obj) {
                B startFuncResult = startFunc.apply(obj);
                return endFunc.apply(startFuncResult);
            }
        };
    }

    /**
     * Join a group of functions with the input and output being the same together
     *
     * @param funcs input functions
     * @param <A>   Input and output object
     * @return a Function that takes A as both input and output
     */
    public static <A> IFunction<A, A> join(final Iterable<? extends IFunction<A, ? extends A>> funcs) {
        return new IFunction<A, A>() {
            @Override
            public A apply(A obj) {
                A result = obj;

                for (IFunction<A, ? extends A> func : funcs) {
                    result = func.apply(result);
                }

                return result;
            }
        };
    }


    /**
     * Join a group of functions with the input and output being the same together
     *
     * @param funcs input functions
     * @param <A>   Input and output object
     * @return a Function that takes A as both input and output
     */
    public static <A> IFunction<A, A> join(final IFunction<A, ? extends A>... funcs) {
        return join(Arrays.asList(funcs));
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


    /**
     * Create a function that does nothing but returning the original object
     *
     * @param <A>   input object
     * @return  original object been input
     */
    public static <A> IFunction<A, A> nothing() {
        return new IFunction<A, A>() {
            @Override
            public A apply(A o) {
                return o;
            }
        };
    }
}

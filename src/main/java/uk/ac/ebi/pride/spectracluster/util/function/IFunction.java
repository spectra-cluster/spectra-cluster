package uk.ac.ebi.pride.spectracluster.util.function;

import java.io.Serializable;

/**
 * Generate an output value for a given input value
 *
 * Implementations of Function interface are expected to be side-effect free and consistent with equals
 * if A == B, then function.apply(A).equals(function.apply(B))
 *
 * Input should be left unchanged
 *
 * @author Rui Wang
 * @version $Id$
 */
public interface IFunction<A, B> extends Serializable{

    B apply(A o);
}

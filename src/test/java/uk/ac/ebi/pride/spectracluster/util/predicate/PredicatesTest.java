package uk.ac.ebi.pride.spectracluster.util.predicate;

import org.junit.Test;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

public class PredicatesTest {

    @Test
    public void testAnd() throws Exception {
        IPredicate<Object> andTrue = Predicates.and(Predicates.alwaysTrue(), Predicates.alwaysTrue(), Predicates.alwaysTrue(), Predicates.alwaysTrue());
        assertTrue(andTrue.apply(new Object()));

        IPredicate<Object> andFalse = Predicates.and(Predicates.alwaysTrue(), Predicates.alwaysTrue(), Predicates.alwaysTrue(), Predicates.alwaysFalse());
        assertFalse(andFalse.apply(new Object()));
    }

    @Test
    public void testOr() throws Exception {
        IPredicate<Object> andTrue = Predicates.or(Predicates.alwaysFalse(), Predicates.alwaysFalse(), Predicates.alwaysFalse(), Predicates.alwaysFalse());
        assertFalse(andTrue.apply(new Object()));

        IPredicate<Object> andFalse = Predicates.or(Predicates.alwaysTrue(), Predicates.alwaysTrue(), Predicates.alwaysTrue(), Predicates.alwaysFalse());
        assertTrue(andFalse.apply(new Object()));
    }

    @Test
    public void testCondition() throws Exception {
        IPredicate<Integer> condition = Predicates.condition(o -> o > 1, o -> ++o);

        assertTrue(condition.apply(1));
        assertFalse(condition.apply(0));
    }

    @Test
    public void testAlwaysTrue() throws Exception {
        IPredicate<Object> predicate = Predicates.alwaysTrue();
        assertTrue(predicate.apply(new Object()));
    }

    @Test
    public void testAlwaysFalse() throws Exception {
        IPredicate<Object> predicate = Predicates.alwaysFalse();
        assertFalse(predicate.apply(new Object()));
    }
}
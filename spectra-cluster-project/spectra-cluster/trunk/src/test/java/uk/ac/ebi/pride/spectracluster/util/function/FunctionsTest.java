package uk.ac.ebi.pride.spectracluster.util.function;

import org.junit.Test;
import uk.ac.ebi.pride.spectracluster.util.predicate.IPredicate;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

public class FunctionsTest {

    @Test
    public void testJoin() throws Exception {
        IFunction<Integer, Integer> minusTwo = Functions.join(new MinusOne(), new MinusOne());

        assertEquals(10, minusTwo.apply(12).intValue());
    }

    @Test
    public void testJoinIterable() throws Exception {
        IFunction<Integer, Integer> minusThree = Functions.join(new MinusOne(), new MinusOne(), new MinusOne());

        assertEquals(10, minusThree.apply(13).intValue());
    }

    @Test
    public void testCondition() throws Exception {
        IFunction<String, Integer> function = new IFunction<String, Integer>() {

            @Override
            public Integer apply(String o) {
                return o.length();
            }
        };

        IPredicate<Integer> predicate = new IPredicate<Integer>() {

            @Override
            public boolean apply(Integer o) {
                return o > 10;
            }
        };

        IFunction<String, Integer> condition = Functions.condition(function, predicate);

        assertEquals(11, condition.apply("KKKKKKKKKKK").intValue());
        assertNull(condition.apply("KKKKKKKKK"));
    }


    private static class MinusOne implements IFunction<Integer, Integer> {

        @Override
        public Integer apply(Integer o) {
            return --o;
        }
    }
}
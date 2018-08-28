package uk.ac.ebi.pride.spectracluster.util;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * com.lordjoe.algorithms.CountedString
 * holds a string and a count - immutable with good hash,equals and compateTO
 * sort by high count then alphabetically
 * User: Steve
 * Date: 7/11/13
 */
public class CountedString implements Comparable<CountedString> {


    /**
     * return an array of CountedString  sorted to by occurance count then
     * alphabetically
     *
     * @param lst !null list of strings - no nulls duplicates OK
     * @return
     */
    public static CountedString[] getCountedStrings(List<String> lst) {
        String[] items = lst.toArray(new String[lst.size()]);
        Arrays.sort(items);
        String current = "";
        int currentCount = 0;
        List<CountedString> holder = new ArrayList<>();

        //noinspection ForLoopReplaceableByForEach
        for (int i = 0; i < items.length; i++) {
            String peptide = items[i];
            if (peptide.equals(current)) {
                currentCount++;
            } else {
                if (!current.isEmpty())    // ignore empty
                    holder.add(new CountedString(current, currentCount));

                currentCount = 1;
                current = peptide;
            }
        }
        if (currentCount > 0) {
            if (!current.isEmpty())    // ignore empty
                holder.add(new CountedString(current, currentCount));

        }
        CountedString[] ret = new CountedString[holder.size()];
        holder.toArray(ret);
        Arrays.sort(ret); // sort by count then alpha
        return ret;
    }


    /**
     * the strings in a list sorted first by occurance count then   alphabetically
     *
     * @param lst !null list of strings - no nulls duplicates OK
     * @return as above
     */
    public static String[] getStringsByOccurance(List<String> lst) {
        CountedString[] sortedByCOubt = getCountedStrings(lst);
        String[] ret = Arrays.stream(sortedByCOubt).map(CountedString::getValue).toArray(String[]::new);
        return ret;
    }


    private final String value;
    private final int count;

    public CountedString(final String pValue, final int pCount) {
        value = pValue;
        count = pCount;
    }

    public String getValue() {
        return value;
    }

    public int getCount() {
        return count;
    }

    /**
     * sort by count then value
     *
     * @param o
     * @return
     */
    @Override
    public int compareTo(final CountedString o) {
        int count = getCount();
        int ocount = o.getCount();
        // High count first
        if (count != ocount)
            return count > ocount ? -1 : 1;
        return getValue().compareTo(o.getValue());
    }

    @Override
    public String toString() {
        return getValue() + ":" + getCount();
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        final CountedString that = (CountedString) o;

        if (count != that.count) return false;
        //noinspection SimplifiableIfStatement,PointlessBooleanExpression,ConstantConditions,RedundantIfStatement
        if (!value.equals(that.value))
            return false;

        return true;
    }

    @Override
    public int hashCode() {
        int result = value.hashCode();
        result = 31 * result + count;
        return result;
    }
}

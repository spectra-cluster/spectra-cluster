package uk.ac.ebi.pride.spectracluster.util.binner;

/**
 * Values may be placed in adjacent bins
 * User: Steve
 * Date: 6/28/13
 */
public interface IWideBinner extends IBinner {

    /**
     * give a list of all bins that the value may be assigned to
     * @param value  value to test
     * @return !null array of bins
     */
    public int[] asBins(double value);


    /**
     * Describe the assigned bins
     * @param value
     * @return either a valid bin number or  null if  isOverflowBinned() is false and the
     *     data is outside the range handled
     */
    public String[] asBinStrings(double value);

}

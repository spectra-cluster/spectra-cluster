package uk.ac.ebi.pride.spectracluster.util.binner;

import java.io.Serializable;

/**
 * Interface describing  a class that divides doubles into bins
 *
 * @author Steve Lewis
 */
public interface IBinner extends Serializable{

    /**
     * place the value into a bin between getMinBin()   and getMaxBin()
     * values outside the range are handled as described below
     *
     * @param value
     * @return either a valid bin number or -1 if isOverflowBinned() is false and the
     * data is outside the range handled
     */
    int asBin(double value);

    /**
     * Describe the assigned bin
     *
     * @param value
     * @return either a valid bin number or  null if  isOverflowBinned() is false and the
     * data is outside the range handled
     */
    String asBinString(double value);

    /**
     * @param bin between
     * @return a number which when sent to asBin will return bin
     * @throws IllegalArgumentException if no such bin is possible
     */
    double fromBin(int bin) throws IllegalArgumentException;


    /**
     * minimum value handed - values below this may be binned as -1 or
     * getMinBin() depending in isOverflowBinned()
     *
     * @return as above
     */
    double getMinValue();

    /**
     * maximum value handed - values below this may be binned as -1 or
     * getMaxBin() depending in isOverflowBinned()
     *
     * @return as above
     */
    double getMaxValue();


    /**
     * minimum bin value - this is almost always 0
     *
     * @return as above
     */
    int getMinBin();

    /**
     * maximum bin value - bins are always LESS than this
     * an array of size getMaxBin() - getMinBin() will hold all legal bins
     *
     * @return as above
     */
    int getMaxBin();

    /**
     * return the total number bins  usually this is the same as getMaxBin
     *
     * @return
     */
    int getNumberBins();


    /**
     * if true values outside getMinValue() .. getMaxValue() are
     * assigned to the highest and lowest bins - otherwise these values return
     * -1
     *
     * @return
     */
    boolean isOverflowBinned();

    /**
     * return this binner but with bins offset by half a bin
     *
     * @return
     */
    IBinner offSetHalf();
}

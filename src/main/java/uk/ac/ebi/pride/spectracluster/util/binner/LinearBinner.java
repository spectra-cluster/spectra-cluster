package uk.ac.ebi.pride.spectracluster.util.binner;

/**
 * Implementation of IBinner as a linear set of bins
 *
 * @author Steve Lewis
 * @date 11/05/13
 */
public class LinearBinner implements IBinner {
    private final double minValue;
    private final double maxValue;
    private final int minBin;
    private final double binSize;
    private final int numberBins;
    private final boolean overFlowBinned;


    // minValue is 0 by default when it is not set

    public LinearBinner(double maxValue, double binSize) {
        this(maxValue, binSize, 0, true, 0);
    }


    public LinearBinner(double maxValue, double binSize, double minValue) {
        this(maxValue, binSize, minValue, true, 0);
    }

    public LinearBinner(double maxValue, double binSize, double minValue, boolean overFlowBinned) {
        this(maxValue, binSize, minValue, overFlowBinned, 0);
    }

    /**
     * create a binner
     *
     * @param maxValue       maximum value
     * @param minValue       minimum value
     * @param binSize        size of bin
     * @param overFlowBinned if true outside range  is binned
     * @param minBin         minimum bin value - usually 0
     */
    public LinearBinner(double maxValue, double binSize, double minValue, boolean overFlowBinned, int minBin) {
        if (maxValue <= minValue)
            throw new IllegalArgumentException("bad bins");
        if (binSize <= 0)
            throw new IllegalArgumentException("bad bins");

        this.minValue = minValue;
        this.maxValue = maxValue;
        this.binSize = binSize;
        this.overFlowBinned = overFlowBinned;
        this.minBin = minBin;


        double del = maxValue - minValue;
        double nb = del / binSize;

        // when rounding a double to an integer add 0.5 to round up, because by default 2.99 turns into 2
        numberBins = (int) (nb + 0.5);
    }

    /**
     * place the value into a bin between getMinBin()   and getMaxBin()
     * values outside the range are handled as described below
     *
     * @param value
     * @return either a valid bin number or -1 if  isOverflowBinned() is false and the
     * data is outside the range handled
     */
    @Override
    public int asBin(double value) {
        double minValue = getMinValue();
        int minBin = getMinBin();
        if (value < minValue) {
            if (isOverflowBinned())
                return minBin;
            else
                return -1; // out of range
        }

        if (value > getMaxValue()) {
            if (isOverflowBinned())
                // -1 is convention for starting with bin 0
                return getMaxBin();
            else
                return -1; // out of range
        }

        double binSize = getBinSize();
        double val = value - minValue;
        int bin = (int) (val / binSize);
        return bin + minBin;
    }

    /**
     * Describe the assigned bin
     *
     * @param value
     * @return either a valid bin number or  null if  isOverflowBinned() is false and the
     * data is outside the range handled
     */
    @Override
    public String asBinString(final double value) {
        int bin = asBin(value);
        if (bin == -1)
            return null;

        return formatBin(bin);

    }

    /**
     * turn a bin into a string
     *
     * @param pBin
     * @return
     */
    protected String formatBin(final int pBin) {
        if (pBin == -1)
            return null;
        double bnv = fromBin(pBin);
        double minv = bnv - getBinSize() / 2;
        double maxv = bnv + getBinSize() / 2;
        return formatBinValue(minv) + "-" + formatBinValue(maxv);
    }

    protected String formatBinValue(double value) {
        return String.format("%10.3f", value).trim();
    }

    public double getBinSize() {
        return binSize;
    }

    /**
     * @param bin between
     * @return a number which when sent to asBin will return bin
     * @throws IllegalArgumentException if no such bin is possible
     */
    @Override
    public double fromBin(int bin) throws IllegalArgumentException {
        if (bin < -1)
            throw new IllegalArgumentException("Illegal bin " + bin);
        if (bin == -1) {
            if (!isOverflowBinned())
                return getMinValue() - 1;
            else
                throw new IllegalArgumentException("Illegal bin " + bin);
        }
        if (bin < getMinBin() || bin > getMaxBin())
            throw new IllegalArgumentException("Illegal bin " + bin);

        // return the bin midpoint
        return getMinValue() + ((bin - getMinBin()) * getBinSize()) + getBinSize() / 2;

    }

    /**
     * minimum value handed - values below this may be binned as -1 or
     * getMinBin() depending in isOverflowBinned()
     *
     * @return as above
     */
    @Override
    public double getMinValue() {
        return minValue;
    }

    /**
     * maximim value handed - values below this may be binned as -1 or
     * getMaxBin() depending in isOverflowBinned()
     *
     * @return as above
     */
    @Override
    public double getMaxValue() {
        return maxValue;
    }

    /**
     * minimum bin value - this is almost always 0
     *
     * @return as above
     */
    @Override
    public int getMinBin() {
        return minBin;
    }

    /**
     * Maximum bin value - bins are always LESS than this
     * an array of size getMaxBin() - getMinBin() will hold all legal bins
     *
     * @return as above
     */
    @Override
    public int getMaxBin() {
        return getMinBin() + getNumberBins() - 1;
    }

    /**
     * if true values outside getMinValue() .. getMaxValue() are
     * assigned to the highest and lowest bins - otherwise these values return
     * -1
     *
     * @return
     */
    @Override
    public boolean isOverflowBinned() {
        return overFlowBinned;
    }

    /**
     * return this binner but with bins offset by half a bin
     *
     * @return
     */
    @Override
    public IBinner offSetHalf() {
        return new LinearBinner(getMaxValue(), getBinSize(), getMinValue() - getBinSize() / 2, isOverflowBinned(), getMinBin());
    }

    /**
     * return the total number bins
     *
     * @return
     */
    public int getNumberBins() {
        return numberBins;
    }
}

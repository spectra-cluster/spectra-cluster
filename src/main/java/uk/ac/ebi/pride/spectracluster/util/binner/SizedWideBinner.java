package uk.ac.ebi.pride.spectracluster.util.binner;

import java.util.Arrays;

/**
 *
 * User: Steve
 * Date: 6/28/13
 */
public class SizedWideBinner extends LinearBinner implements IWideBinner {

    private final double overlapWidth;

    public SizedWideBinner(final double maxValue, final double binSize, final double minValue, final double overlapWidth) {
        this(maxValue, binSize, minValue, overlapWidth, false);
    }

    public SizedWideBinner(final double maxValue, final double binSize, final double minValue, final double overlapWidth, final boolean overFlowBinned) {
        super(maxValue, binSize, minValue, overFlowBinned);
        this.overlapWidth = overlapWidth;
        if (overlapWidth > binSize)
            throw new IllegalArgumentException("overlapWidth must be less than binSize");
    }

    protected double getOverlapWidth() {
        return overlapWidth;
    }

    /**
     * give a list of all bins that the value may be assigned to
     *
     * @param value value to test
     * @return !null array of bins
     */
    @Override
    public int[] asBins(final double value) {
        int mainBin = asBin(value);

        double overlapWidth = getOverlapWidth();
        double lowValue = value - overlapWidth;
        lowValue = Math.max(lowValue, getMinValue());
        int lowBin = asBin(lowValue);
        if (lowBin < mainBin) {
            return new int[]{lowBin, mainBin};
        }

        double highValue = value + overlapWidth;
        highValue = Math.min(highValue, getMaxValue());
        int highBin = asBin(highValue);
        if (highBin > mainBin) {
            return new int[]{mainBin, highBin};
        }

        return new int[]{mainBin};
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
        double v = getBinSize() / 2 + getOverlapWidth();
        double minv = bnv - v;
        double maxv = bnv + v;
        return formatBinValue(minv) + "-" + formatBinValue(maxv);
    }


    /**
     * Describe the assigned bins
     *
     * @param value
     * @return either a valid bin number or  null if  isOverflowBinned() is false and the
     * data is outside the range handled
     */
    @Override
    public String[] asBinStrings(final double value) {
        int[] bins = asBins(value);
        return Arrays.stream(bins)
                .mapToObj(this::formatBin)
                .toArray(String[]::new);
    }

    /**
     * return this binner but with bins offset by half a bin
     *
     * @return
     */
    @Override
    public IBinner offSetHalf() {
        return new SizedWideBinner(getMaxValue(), getBinSize(), getMinValue() - getBinSize() / 2, getOverlapWidth(), isOverflowBinned());

    }

}

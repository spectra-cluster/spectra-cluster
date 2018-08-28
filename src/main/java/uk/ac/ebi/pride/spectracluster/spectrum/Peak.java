package uk.ac.ebi.pride.spectracluster.spectrum;

import uk.ac.ebi.pride.spectracluster.util.MZIntensityUtilities;

/**
 * @author Steve Lewis
 * @author Rui Wang
 * @version $Id$
 */
public class Peak implements IPeak {

    private final float massChargeRatio;
    private final float intensity;
    private final int count;

    public Peak(float massChargeRatio, float intensity) {
        this(massChargeRatio, intensity, 1);
    }

    public Peak(float massChargeRatio, float intensity, int count) {
        this.massChargeRatio = massChargeRatio;
        this.intensity = intensity;
        this.count = count;
    }

    /**
     * copy constructor
     *
     * @param copied
     */
    public Peak(IPeak copied) {
        this.massChargeRatio = copied.getMz();
        this.intensity = copied.getIntensity();
        this.count = copied.getCount();
    }

    public float getMz() {
        return massChargeRatio;
    }

    public float getIntensity() {
        return intensity;
    }

    public int getCount() {
        return count;
    }

    @Override
    public int compareTo(IPeak o) {
        float mz = getMz();
        float omz = o.getMz();

        int ret =  Float.compare(mz, omz);
        if(ret != 0) return ret;

        float intent = getIntensity();
        float ointent = o.getIntensity();

        ret = Float.compare(intent, ointent);
        if(ret != 0) return ret;
        // test count or return 0
        return 0; // same

    }

//    This compareTo method can cause exceptions in sorting in Java 7, due to a change in sorting related algorithm
//    @Override
//    public int compareTo(IPeak o) {
//        if (Math.abs(getMz() - o.getMz()) > MZIntensityUtilities.SMALL_MZ_DIFFERENCE)
//            return CompareTo.compare(getMz(), o.getMz());
//        if (Math.abs(getIntensity() - o.getIntensity()) > MZIntensityUtilities.SMALL_INTENSITY_DIFFERENCE)
//            return Double.compare(getIntensity(), o.getIntensity());
//        return 0;
//    }

    /**
     * like equals but weaker - says other is equivalent to this
     *
     * @param other possible null other object
     * @return true if other is "similar enough to this"
     */
    @Override
    public boolean equivalent(IPeak other) {
        if (getCount() != other.getCount())
            return false;

        Float mzDiff = Math.abs(other.getMz() - getMz());
        if (mzDiff > MZIntensityUtilities.SMALL_MZ_DIFFERENCE)
            return false;

        Float intensDiff = Math.abs(other.getIntensity() - getIntensity());
        return !(intensDiff > 0.001);
    }

    /**
     * return exactly what an MGF would use
     *
     * @return
     */
    @Override
    public String toString() {
        String mz = String.format("%10.5f", getMz()).trim();
        String intensity = String.format("%8.2f", getIntensity()).trim();
        return "m/z = " + mz + ", intensity = " + intensity + ", count = " + getCount();
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        final Peak peak = (Peak) o;

        return Float.compare(peak.intensity, intensity) == 0 && Float.compare(peak.massChargeRatio, massChargeRatio) == 0;

    }

    @Override
    public int hashCode() {
        int result = (massChargeRatio != +0.0f ? Float.floatToIntBits(massChargeRatio) : 0);
        result = 31 * result + (intensity != +0.0f ? Float.floatToIntBits(intensity) : 0);
        return result;
    }
}

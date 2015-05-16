package uk.ac.ebi.pride.spectracluster.spectrum;

import org.junit.Assert;
import org.junit.Test;

/**
 * Created by jg on 13.05.15.
 */
public class TestMasses {
    @Test
    public void testMasses() {
        float mz = 100.123F;
        int charge = 3;

        float mass = Masses.getMonoisotopicMass(mz, charge);
        float calcMz = Masses.getMz(mass, charge);
        Assert.assertEquals(calcMz, mz, 0.00001);
    }
}

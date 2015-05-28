package uk.ac.ebi.pride.spectracluster.util.predicate.spectrum_comparison;

import junit.framework.Assert;
import org.junit.Before;
import org.junit.Test;
import uk.ac.ebi.pride.spectracluster.io.ParserUtilities;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Created by jg on 28.05.15.
 */
public class ShareMajorPeaksPredicateTest {
    private List<ISpectrum> testSpectra;

    @Before
    public void setUp() throws Exception {
        File testFile = new File(ShareMajorPeaksPredicateTest.class.getClassLoader().getResource("spectra_400.0_4.0.mgf").toURI());
        testSpectra = new ArrayList<ISpectrum>();
        ISpectrum[] readSpectra = ParserUtilities.readMGFScans(testFile);

        for (ISpectrum s : readSpectra)
            testSpectra.add(s);

    }

    @Test
    public void testSharedPeaks() throws Exception {
        ShareMajorPeaksPredicate predicate = new ShareMajorPeaksPredicate(5);

        Assert.assertFalse(predicate.apply(testSpectra.get(0), testSpectra.get(1)));
        Assert.assertTrue(predicate.apply(testSpectra.get(0), testSpectra.get(3)));
        Assert.assertTrue(predicate.apply(testSpectra.get(0), testSpectra.get(5)));
        Assert.assertTrue(predicate.apply(testSpectra.get(0), testSpectra.get(7)));
        Assert.assertFalse(predicate.apply(testSpectra.get(0), testSpectra.get(2)));
    }
}

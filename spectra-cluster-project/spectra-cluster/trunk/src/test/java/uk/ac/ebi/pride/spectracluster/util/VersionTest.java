package uk.ac.ebi.pride.spectracluster.util;

import org.junit.*;

/**
 * uk.ac.ebi.pride.spectracluster.util.VersionTest
 *
 * @author Steve Lewis
 * @date 05/06/2014
 */
public class VersionTest {
    // Test that version is NOT the default value
    @Test
    public void testVersion()
    {
        final String version = Version.version;
        Assert.assertFalse("This only works in a released and tagged version",Version.SNAPSHOT.equals(version));
    }
}

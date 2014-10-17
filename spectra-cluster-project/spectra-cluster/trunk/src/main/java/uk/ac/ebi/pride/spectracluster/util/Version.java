package uk.ac.ebi.pride.spectracluster.util;

import java.io.File;
import java.io.FileInputStream;
import java.net.URL;
import java.util.Properties;

/**
 * uk.ac.ebi.pride.spectracluster.util.Version
 *
 * @author Steve Lewis
 * @author Rui Wang
 * @date 05/06/2014
 */
public class Version {

    public static final String SNAPSHOT = "1.0.1-SNAPSHOT";
    public static String version = SNAPSHOT;

    static {
        // get the location of version property file
        URL versionPropertyFileUrl = Version.class.getClassLoader().getResource("version.properties");

        // load version property file
        final Properties versionProperties = new Properties();
        try {
            versionProperties.load(new FileInputStream(new File(versionPropertyFileUrl.toURI())));
            final String property = versionProperties.getProperty("algorithm.version");
            if (!property.startsWith("$")) {
                version = property;
            }
        } catch (Exception e) {
            version = SNAPSHOT;
        }

    }


}

package uk.ac.ebi.pride.spectracluster.cdf;

import uk.ac.ebi.pride.spectracluster.similarity.CombinedFisherIntensityTest;
import uk.ac.ebi.pride.spectracluster.similarity.FrankEtAlDotProduct;
import uk.ac.ebi.pride.spectracluster.similarity.ISimilarityChecker;
import uk.ac.ebi.pride.spectracluster.util.Defaults;

import java.io.BufferedReader;
import java.io.InputStreamReader;

/**
 * Created by jg on 05.05.15.
 */
public class CumulativeDistributionFunctionFactory {
    private CumulativeDistributionFunctionFactory() {

    }

    /**
     * Returns the default cumulative distribution function (CDF). If the CDF is
     * set in the Defaults class, this set CDF is returned. Otherwise the matching
     * CDF is loaded from the resource file.
     * @param similatiryClass
     * @return
     * @throws Exception
     */
    public static CumulativeDistributionFunction getDefaultCumlativeDistributionFunctionForSimilarityMetric(Class similatiryClass) throws Exception {
        if (Defaults.getCumulativeDistributionFunction() != null) {
            return Defaults.getCumulativeDistributionFunction();
        }
        else {
            return getCumulativeDistributionFunctionForSimilarityMetric(similatiryClass);
        }
    }

    /**
     * Returns the matching cumulative distribution function from the matching
     * resource file. If no resource file exists for the passed similarity
     * checker class an Exception is thrown.
     * @param similarityCheckerClass
     * @return
     * @throws Exception Thrown if no CDF resource exists for the passed similarity checker class.
     */
    public static CumulativeDistributionFunction getCumulativeDistributionFunctionForSimilarityMetric(Class similarityCheckerClass) throws Exception {
        if (similarityCheckerClass == CombinedFisherIntensityTest.class) {
            return getCumulativeDistributionFunctionForResource("cumulative.cdf.tsv");
        }
        if (similarityCheckerClass == FrankEtAlDotProduct.class) {
            return getCumulativeDistributionFunctionForResource("dot.cdf.tsv");
        }

        throw new Exception("No cumulative distribution function defined for " + similarityCheckerClass);
    }

    /**
     * Reads the cumulative distribution function from a resource and builds the
     * matching CumulativeDistributionFunction class from it.
     * @param resource
     * @return
     * @throws Exception
     */
    private static CumulativeDistributionFunction getCumulativeDistributionFunctionForResource(String resource) throws Exception {
        BufferedReader reader = new BufferedReader(new InputStreamReader(CumulativeDistributionFunctionFactory.class.getClassLoader().getResourceAsStream(resource)));

        StringBuilder definitionString = new StringBuilder();

        String line;

        while ((line = reader.readLine()) != null)
            definitionString.append(line).append("\n");

        return CumulativeDistributionFunction.fromString(definitionString.toString());
    }
}

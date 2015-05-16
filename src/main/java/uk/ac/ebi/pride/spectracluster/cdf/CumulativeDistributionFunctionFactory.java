package uk.ac.ebi.pride.spectracluster.cdf;

import uk.ac.ebi.pride.spectracluster.similarity.CombinedFisherIntensityTest;
import uk.ac.ebi.pride.spectracluster.similarity.ISimilarityChecker;

import java.io.BufferedReader;
import java.io.InputStreamReader;

/**
 * Created by jg on 05.05.15.
 */
public class CumulativeDistributionFunctionFactory {
    private CumulativeDistributionFunctionFactory() {

    }

    public static CumulativeDistributionFunction getCumulativeDistributionFunctionForSimilarityMetric(Class similarityCheckerClass) throws Exception {
        if (similarityCheckerClass == CombinedFisherIntensityTest.class) {
            return getCumulativeDistributionFunctionForResource("cumulative.cdf.tsv");
        }

        throw new Exception("No cumulative distribution function defined for " + similarityCheckerClass);
    }

    private static CumulativeDistributionFunction getCumulativeDistributionFunctionForResource(String resource) throws Exception {
        BufferedReader reader = new BufferedReader(new InputStreamReader(CumulativeDistributionFunctionFactory.class.getClassLoader().getResourceAsStream(resource)));

        StringBuilder definitionString = new StringBuilder();

        String line;

        while ((line = reader.readLine()) != null)
            definitionString.append(line).append("\n");

        return CumulativeDistributionFunction.fromString(definitionString.toString());
    }
}

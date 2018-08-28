package uk.ac.ebi.pride.spectracluster.similarity;

import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.stat.correlation.KendallsCorrelation;

import java.io.Serializable;

/**
 * This code is licensed under the Apache License, Version 2.0 (the
 * "License"); you may not use this file except in compliance
 * with the License.  You may obtain a copy of the License at
 * <p>
 * http://www.apache.org/licenses/LICENSE-2.0
 * <p>
 * ==Overview==
 * <p>
 * This class
 * <p>
 * Created by ypriverol (ypriverol@gmail.com) on 07/11/2017.
 */
public class SerializableKendallsCorrelation extends KendallsCorrelation implements Serializable{

    public SerializableKendallsCorrelation(double[][] data) {
        super(data);
    }

    public SerializableKendallsCorrelation(RealMatrix matrix) {
        super(matrix);
    }

    public SerializableKendallsCorrelation() {

    }
}

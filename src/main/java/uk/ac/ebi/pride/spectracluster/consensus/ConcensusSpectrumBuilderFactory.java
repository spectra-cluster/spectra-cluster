package uk.ac.ebi.pride.spectracluster.consensus;

import java.io.Serializable;

/**
 * uk.ac.ebi.pride.spectracluster.consensus.ConcensusSpectrumBuilderFactory
 * User: Steve
 * Date: 8/7/13
 */
public interface ConcensusSpectrumBuilderFactory extends Serializable{

    /**
     * build a new instance of the cpectrum builder
     *
     * @return !null instance
     */
    IConsensusSpectrumBuilder getConsensusSpectrumBuilder();

}

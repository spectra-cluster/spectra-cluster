package uk.ac.ebi.pride.spectracluster.util;

import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;

import java.util.List;

/**
 * uk.ac.ebi.pride.spectracluster.spectrum.ConcensusSpectraItems
 *
 * @author Steve Lewis
 */
public class ConsensusSpectraItems {

    private ISpectrum concensus;
    private List<ISpectrum> spectra;

    public ISpectrum getConcensus() {
        return concensus;
    }

    public void setConcensus(ISpectrum concensus) {
        this.concensus = concensus;
    }

    public List<ISpectrum> getSpectra() {
        return spectra;
    }

    public void setSpectra(List<ISpectrum> spectra) {
        this.spectra = spectra;
    }
}

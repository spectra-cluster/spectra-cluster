package uk.ac.ebi.pride.spectracluster.io;

import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;

/**
 * Appender interface for Spectrum
 *
 * @author Rui Wang
 * @version $Id$
 */
public interface ISpectrumAppender {

    /**
     * Add spectrum
     *
     * @param out       !!null open appendable
     * @param spectrum  !null cluster
     * @param otherData any other data - implementation specific and usually blank
     */
    void appendSpectrum(Appendable out, ISpectrum spectrum, Object... otherData) throws AppenderException;

    /**
     * add whatever happens at the start
     *
     * @param out       !null open appendable
     * @param otherData any other data - implementation specific and usually blank
     */
    void appendStart(Appendable out, Object... otherData) throws AppenderException;

    /**
     * add whatever happens at the end
     *
     * @param out       !null open appendable
     * @param otherData any other data - implementation specific and usually blank
     */
    void appendEnd(Appendable out, Object... otherData) throws AppenderException;
}

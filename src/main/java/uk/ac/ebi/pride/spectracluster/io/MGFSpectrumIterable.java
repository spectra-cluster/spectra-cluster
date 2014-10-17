package uk.ac.ebi.pride.spectracluster.io;

import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;

import java.io.*;
import java.util.Iterator;

/**
 * uk.ac.ebi.pride.spectracluster.spectrum.MGFSpectrumIterable
 * class to read spectra from a MGF file one at a time
 * This allows clusters to be processed from a large file without
 * blowing memory limits* User: Steve
 * Date: 7/19/13
 */
public class MGFSpectrumIterable implements Iterable<ISpectrum> {

    protected static LineNumberReader fileToLineNumberReader(File f) {
        try {
            return new LineNumberReader(new FileReader(f));
        } catch (FileNotFoundException e) {
            throw new RuntimeException(e);

        }

    }

    private final LineNumberReader reader;
    private final MGFSpectrumIterator one_time_iterator;
    private ISpectrum nextSpectrum;

    /**
     * build with an existing readable file
     *
     * @param f !null existing non-directory mgf file
     */
    public MGFSpectrumIterable(File f) {
        this(fileToLineNumberReader(f));
    }

    /**
     * build with input stream
     *
     * @param rdr !null open input stream
     */
    public MGFSpectrumIterable(InputStream rdr) {
        this(new LineNumberReader(new InputStreamReader(rdr)));
    }

    /**
     * build with reader
     *
     * @param rdr !null open reader
     */
    public MGFSpectrumIterable(Reader rdr) {
        this(new LineNumberReader(rdr));
    }


    /**
     * build with LineNumberReader
     *
     * @param rdr !null open LineNumberReader
     */
    public MGFSpectrumIterable(LineNumberReader rdr) {
        reader = new LineNumberReader(rdr);
        nextSpectrum = ParserUtilities.readMGFScan(reader);
        one_time_iterator = new MGFSpectrumIterator();

    }


    /**
     * Returns an iterator over a set of elements of type T.
     *
     * @return an Iterator.
     */
    @Override
    public Iterator<ISpectrum> iterator() {
        return one_time_iterator;
    }


    protected class MGFSpectrumIterator implements Iterator<ISpectrum> {
        /**
         * Returns <tt>true</tt> if the iteration has more elements. (In other
         * words, returns <tt>true</tt> if <tt>next</tt> would return an element
         * rather than throwing an exception.)
         *
         * @return <tt>true</tt> if the iterator has more elements.
         */
        @Override
        public boolean hasNext() {
            return nextSpectrum != null;
        }

        /**
         * Returns the next element in the iteration.
         *
         * @return the next element in the iteration.
         * @throws java.util.NoSuchElementException iteration has no more elements.
         */
        @Override
        public ISpectrum next() {
            ISpectrum ret = nextSpectrum;
            nextSpectrum = ParserUtilities.readMGFScan(reader);
            return ret;
        }

        @Override
        public void remove() {
            throw new UnsupportedOperationException("Not allowed");
        }
    }


}

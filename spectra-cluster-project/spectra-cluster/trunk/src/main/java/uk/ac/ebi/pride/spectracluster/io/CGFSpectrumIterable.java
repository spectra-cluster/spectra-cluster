package uk.ac.ebi.pride.spectracluster.io;

import uk.ac.ebi.pride.spectracluster.cluster.ICluster;

import java.io.*;
import java.util.Iterator;

/**
 * uk.ac.ebi.pride.spectracluster.spectrum.CGFSpectrumIterable
 * class to read clusters from a CGF file one at a time
 * This allows clusters to be processed from a large file without
 * blowing memory limits
 * User: Steve
 * Date: 7/19/13
 */
public class CGFSpectrumIterable implements Iterable<ICluster> {

    protected static LineNumberReader fileToLineNumberReader(File f) {
        try {
            return new LineNumberReader(new FileReader(f));
        } catch (FileNotFoundException e) {
            throw new RuntimeException(e);

        }

    }

    private final LineNumberReader reader;
    private final MGFSpectrumIterator one_time_iterator;
    private ICluster nextCluster;

    /**
     * build with an existing readable file
     *
     * @param f !null existing non-directory mgf file
     */
    @SuppressWarnings("UnusedDeclaration")
    public CGFSpectrumIterable(File f) {
        this(fileToLineNumberReader(f));
    }

    /**
     * build with input stream
     *
     * @param rdr !null open input stream
     */
    @SuppressWarnings("UnusedDeclaration")
    public CGFSpectrumIterable(InputStream rdr) {
        this(new LineNumberReader(new InputStreamReader(rdr)));
    }

    /**
     * build with reader
     *
     * @param rdr !null open reader
     */
    @SuppressWarnings("UnusedDeclaration")
    public CGFSpectrumIterable(Reader rdr) {
        this(new LineNumberReader(rdr));
    }


    /**
     * build with LineNumberReader
     *
     * @param rdr !null open LineNumberReader
     */
    public CGFSpectrumIterable(LineNumberReader rdr) {
        reader = new LineNumberReader(rdr);
        nextCluster = ParserUtilities.readSpectralCluster(reader, null);
        one_time_iterator = new MGFSpectrumIterator();
    }


    /**
     * Returns an iterator over a set of elements of type T.
     *
     * @return an Iterator.
     */
    @Override
    public Iterator<ICluster> iterator() {
        return one_time_iterator;
    }


    protected class MGFSpectrumIterator implements Iterator<ICluster> {
        /**
         * Returns <tt>true</tt> if the iteration has more elements. (In other
         * words, returns <tt>true</tt> if <tt>next</tt> would return an element
         * rather than throwing an exception.)
         *
         * @return <tt>true</tt> if the iterator has more elements.
         */
        @Override
        public boolean hasNext() {
            return nextCluster != null;
        }

        /**
         * Returns the next element in the iteration.
         *
         * @return the next element in the iteration.
         * @throws java.util.NoSuchElementException iteration has no more elements.
         */
        @Override
        public ICluster next() {
            ICluster ret = nextCluster;
            nextCluster = ParserUtilities.readSpectralCluster(reader, null);
            return ret;
        }

        @Override
        public void remove() {
            throw new UnsupportedOperationException("Not allowed");
        }
    }


}

package uk.ac.ebi.pride.spectracluster.io;

/**
 * Runtime exception for appenders
 *
 * @author Rui Wang
 * @version $Id$
 */
public class AppenderException extends RuntimeException {

    public AppenderException() {
    }

    public AppenderException(String message) {
        super(message);
    }

    public AppenderException(String message, Throwable cause) {
        super(message, cause);
    }

    public AppenderException(Throwable cause) {
        super(cause);
    }
}

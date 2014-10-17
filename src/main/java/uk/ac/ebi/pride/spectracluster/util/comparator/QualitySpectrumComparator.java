package uk.ac.ebi.pride.spectracluster.util.comparator;

import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.util.CompareTo;

import java.util.Comparator;

/**
 * @author Rui Wang
 * @version $Id$
 */
public class QualitySpectrumComparator implements Comparator<ISpectrum> {

    /**
     * sort with the highest quality spectra first
     *
     * @param cluster1
     * @param cluster2
     * @return
     */
    @Override
    public int compare(ISpectrum cluster1, ISpectrum cluster2) {
        if (cluster1 == cluster2)
            return 0;

        double q1 = cluster1.getQualityScore();
        double q2 = cluster2.getQualityScore();
        return CompareTo.compare(q2, q1);

    }


}

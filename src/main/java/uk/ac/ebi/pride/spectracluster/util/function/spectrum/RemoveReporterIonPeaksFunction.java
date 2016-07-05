package uk.ac.ebi.pride.spectracluster.util.function.spectrum;

import uk.ac.ebi.pride.spectracluster.spectrum.IPeak;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.spectrum.Masses;
import uk.ac.ebi.pride.spectracluster.spectrum.Spectrum;
import uk.ac.ebi.pride.spectracluster.util.function.IFunction;

import java.util.ArrayList;
import java.util.List;

/**
 * Removes reporter ions associated peaks from the passed spectrum.
 * Currently, the class simply removes the highest number of peaks
 * known for a given quantitation method (ie. iTRAQ 8 peaks for iTRAQ
 * and TMT 10 plex peaks for all TMT experiments).
 *
 * Created by jg on 13.05.15.
 */
public class RemoveReporterIonPeaksFunction implements IFunction<ISpectrum, ISpectrum> {
    private final float fragmentIonTolerance;
    private final REPORTER_TYPE reporterType;

    public enum REPORTER_TYPE {
        ITRAQ,
        TMT,
        ALL
    }

    public RemoveReporterIonPeaksFunction(float fragmentIonTolerance, REPORTER_TYPE reporterType) {
        this.fragmentIonTolerance = fragmentIonTolerance;
        this.reporterType = reporterType;
    }

    /**
     * By default all reporter ions within the defined fragment
     * tolerance are removed.
     * @param fragmentIonTolerance
     */
    public RemoveReporterIonPeaksFunction(float fragmentIonTolerance) {
        this(fragmentIonTolerance, REPORTER_TYPE.ALL);
    }

    @Override
    public ISpectrum apply(ISpectrum o) {

        // get the m/z values
        float[] reporterMzValues = getReporterMz(reporterType);

        // filter the peak list
        List<IPeak> filteredPeakList = new ArrayList<IPeak>();

        for (IPeak peak : o.getPeaks()) {
            final float peakMz = peak.getMz();

            // ignore any peak that could be a neutral loss
            boolean isReporterPeak = false;
            for (float reporterMz : reporterMzValues) {
                if (isWithinRange(reporterMz - fragmentIonTolerance, reporterMz + fragmentIonTolerance, peakMz)) {
                    isReporterPeak = true;
                    break;
                }
            }
            if (!isReporterPeak) {
                filteredPeakList.add(peak);
            }
        }

        ISpectrum filteredSpectrum = new Spectrum(o, filteredPeakList, true);

        return filteredSpectrum;
    }

    /**
     * Returns the m/z values for the reporter ions relevant
     * to the passed reporterType.
     * @param reporterType
     * @return
     */
    private float[] getReporterMz(REPORTER_TYPE reporterType) {
        switch (reporterType) {
            case ITRAQ:
                // 305.1 for iTRAQ 9 was explicitly not added yet
                return new float[] {113.1F, 114.1F, 115.1F, 116.1F, 117.1F, 118.1F, 119.1F, 121.1F};
            case TMT:
                // 230.1697 represents the complete TMT tag
                return new float[] {126.127725F, 127.12476F, 127.131079F, 128.128114F, 128.134433F, 129.131468F,
                                    129.137787F, 130.134822F, 130.141141F, 131.138176F, 230.1697f};
            case ALL:
            default:
                // merge all known reporters
                List<Float> reporterMz = new ArrayList<Float>();
                for (REPORTER_TYPE rt : REPORTER_TYPE.values()) {
                    if (rt == REPORTER_TYPE.ALL)
                        continue;
                    float[] curRtsMz = getReporterMz(rt);
                    for (float mz : curRtsMz)
                        reporterMz.add(mz);
                }

                // change to float[]
                float[] returnVal = new float[reporterMz.size()];
                for (int i = 0; i < reporterMz.size(); i++)
                    returnVal[i] = reporterMz.get(i).floatValue();

                return returnVal;
        }
    }

    private boolean isWithinRange(float min, float max, float value) {
        return (value >= min && value <= max);
    }
}

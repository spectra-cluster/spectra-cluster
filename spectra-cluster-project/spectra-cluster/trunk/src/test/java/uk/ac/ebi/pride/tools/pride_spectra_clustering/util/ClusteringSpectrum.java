package uk.ac.ebi.pride.tools.pride_spectra_clustering.util;

import uk.ac.ebi.pride.tools.jmzreader.model.Spectrum;
import uk.ac.ebi.pride.tools.jmzreader.model.impl.ParamGroup;
import uk.ac.ebi.pride.tools.jmzreader.model.impl.SpectrumImplementation;
import uk.ac.ebi.pride.tools.mgf_parser.model.Ms2Query;

import java.util.*;

/**
 * NOTE extend WatchedClass to look for possible memory leaks
 */
public class ClusteringSpectrum /* extends WatchedClass */ {
    private String id;
    private String peptide;
    private Double precursorMZ;
    private Double precursorIntensity;
    private Integer precursorCharge;
    private List<Peak> peaklist;
    private ParamGroup additional;
    private Integer msLevel;

    public ClusteringSpectrum(Spectrum s) {

        id = s.getId();
        if (s instanceof Ms2Query) {
            id = ((Ms2Query) s).getTitle();
            int index = id.indexOf(",sequence=");
            if (index > -1) {
                peptide = id.substring(index + ",sequence=".length()).trim();
                id = id.substring(0, index).replace("id=", "");
            }
        }

        precursorMZ = s.getPrecursorMZ();
        precursorIntensity = s.getPrecursorIntensity();
        precursorCharge = s.getPrecursorCharge();
        additional = s.getAdditional();
        msLevel = s.getMsLevel();

        peaklist = new ArrayList<Peak>(s.getPeakList().size());

        for (Double mz : s.getPeakList().keySet())
            peaklist.add(new Peak(mz, s.getPeakList().get(mz)));

        Collections.sort(peaklist, PeakIntensityComparator.getInstance());
    }

    public ClusteringSpectrum(String id, Double precursorMZ,
                              Double precursorIntensity, Integer precursorCharge,
                              Map<Double, Double> peaklist, ParamGroup additional, Integer msLevel) {
        this.id = id;
        this.precursorMZ = precursorMZ;
        this.precursorIntensity = precursorIntensity;
        this.precursorCharge = precursorCharge;
        this.additional = additional;
        this.msLevel = msLevel;

        this.peaklist = new ArrayList<Peak>(peaklist.size());

        for (Double mz : peaklist.keySet())
            this.peaklist.add(new Peak(mz, peaklist.get(mz)));

        Collections.sort(this.peaklist, PeakIntensityComparator.getInstance());
    }

    /**
     * Added SLewis because we hand a List of parks not a Map
     *
     * @param id
     * @param precursorMZ
     * @param precursorIntensity
     * @param precursorCharge
     * @param peaklist
     */
    public ClusteringSpectrum(String id,
                              double precursorMZ,
                              double precursorIntensity,
                              int precursorCharge,
                              List<Peak> peaklist) {
        this.id = id;
        this.precursorMZ = precursorMZ;
        this.precursorIntensity = precursorIntensity;
        this.precursorCharge = precursorCharge;
        this.additional = null;
        this.msLevel = 2;

        this.peaklist = new ArrayList<Peak>(peaklist);

        Collections.sort(this.peaklist, PeakIntensityComparator.getInstance());
    }


    public String getPeptide() {
        return peptide;
    }

    public void setPeptide(final String pPeptide) {
        peptide = pPeptide;
    }

    public String getId() {
        return id;
    }

    public Double getPrecursorMZ() {
        return precursorMZ;
    }

    public Double getPrecursorIntensity() {
        return precursorIntensity;
    }

    public Integer getPrecursorCharge() {
        return precursorCharge;
    }

    public List<Peak> getPeaklist() {
        return peaklist;
    }

    public ParamGroup getAdditional() {
        return additional;
    }

    public Spectrum toSpectrum() {
        Map<Double, Double> mapPeaks = new HashMap<Double, Double>(peaklist.size());
        for (Peak p : peaklist)
            mapPeaks.put(p.getMz(), p.getIntensity());
        return new SpectrumImplementation(id, precursorCharge, precursorMZ, precursorIntensity, mapPeaks, msLevel);
    }

    @Override
    public String toString() {
        String ret = id;
        if (peptide != null)
            ret += "," + peptide;
        return ret;
    }
}

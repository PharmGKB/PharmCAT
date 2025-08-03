package org.pharmgkb.pharmcat.definition.model;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.SortedSet;
import com.google.common.base.Preconditions;
import com.google.common.base.Splitter;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;
import org.apache.commons.lang3.StringUtils;
import org.jspecify.annotations.Nullable;
import org.pharmgkb.common.comparator.ChromosomeNameComparator;


/**
 * All the information used to describe a particular location of a variant.
 *
 * @author Ryan Whaley
 */
public class VariantLocus implements Comparable<VariantLocus> {
  public static final Splitter HGVS_NAME_SPLITTER = Splitter.on(";").trimResults();
  @Expose
  @SerializedName("chromosome")
  private final String m_chromosome;
  /**
   * Position as mapped for VCF.
   */
  @Expose
  @SerializedName("position")
  private long m_position;
  /**
   * Position as defined in CPIC.
   */
  @Expose
  @SerializedName("cpicPosition")
  private final long m_cpicPosition;
  @Expose
  @SerializedName("rsid")
  private @Nullable String m_rsid;
  @Expose
  @SerializedName("chromosomeHgvsName")
  private final String m_chromosomeHgvsName;
  private List<String> m_chromosomeHgvsNameList;
  /**
   * Alleles as defined in CPIC.
   */
  @Expose
  @SerializedName("cpicAlleles")
  private SortedSet<String> m_cpicAlleles;
  /**
   * Map of CPIC to VCF alleles.
   */
  @Expose
  @SerializedName("cpicToVcfAlleleMap")
  private Map<String, String> m_cpicToVcfAlleleMap;

  /**
   * VCF ref allele.
   */
  @Expose
  @SerializedName("ref")
  private String m_ref;
  /**
   * VCF alt alleles.
   */
  @Expose
  @SerializedName("alts")
  private List<String> m_alts;


  public VariantLocus(String chromosome, long position, String chromosomeHgvsName) {
    Preconditions.checkNotNull(chromosome);
    Preconditions.checkNotNull(chromosomeHgvsName);
    m_chromosome = chromosome;
    m_position = position;
    m_cpicPosition = position;
    m_chromosomeHgvsName = chromosomeHgvsName;
    m_chromosomeHgvsNameList = HGVS_NAME_SPLITTER.splitToList(m_chromosomeHgvsName);
  }


  /**
   * Gets the name of the chromosome (e.g. "chr12").
   */
  public String getChromosome() {
    return m_chromosome;
  }


  /**
   * Gets the chromosome and VCF position for this variant (e.g. "chr12:234948").
   */
  public String getVcfChrPosition() {
    return m_chromosome + ":" + m_position;
  }


  /**
   * Gets the (start) position on the chromosomal sequence as mapped for VCF.
   */
  public long getPosition() {
    return m_position;
  }

  public void setPosition(long position) {
    m_position = position;
  }


  /**
   * Gets the (start) position on the chromosomal sequence as defined by CPIC.
   */
  public long getCpicPosition() {
    return m_cpicPosition;
  }


  /**
   * The HGVS name used for this location on the chromosomal sequence (NC_).
   */
  public String getChromosomeHgvsName() {
    return m_chromosomeHgvsName;
  }


  /**
   * The identifier to use for this location from dbSNP
   */
  public @Nullable String getRsid() {
    return m_rsid;
  }

  public void setRsid(@Nullable String rsid) {
    m_rsid = rsid;
  }


  /**
   * Gets the alleles as defined in CPIC.
   */
  public SortedSet<String> getCpicAlleles() {
    return m_cpicAlleles;
  }

  public void setCpicAlleles(SortedSet<String> cpicAlleles) {
    m_cpicAlleles = cpicAlleles;
  }


  /**
   * Gets the map of CPIC to VCF alleles.
   */
  public Map<String, String> getCpicToVcfAlleleMap() {
    return m_cpicToVcfAlleleMap;
  }

  public void setCpicToVcfAlleleMap(Map<String, String> cpicToVcfAlleleMap) {
    m_cpicToVcfAlleleMap = cpicToVcfAlleleMap;
  }

  /**
   * Checks if specified {@code allele} is a ref or alt allele.
   */
  public boolean hasVcfAllele(String allele) {
    return m_cpicToVcfAlleleMap.containsValue(allele);
  }


  /**
   * Gets the VCF ref allele.
   */
  public String getRef() {
    return m_ref;
  }

  public void setRef(String ref) {
    Preconditions.checkNotNull(ref);
    m_ref = ref;
  }

  /**
   * Gets the VCF alt alleles.
   */
  public List<String> getAlts() {
    return m_alts;
  }

  public void setAlts(List<String> alts) {
    m_alts = alts;
  }

  public void addAlt(String alt) {
    if (m_alts == null) {
      m_alts = new ArrayList<>();
    }
    m_alts.add(alt);
  }


  public String getHgvsForVcfAllele(String vcfAllele) {
    Preconditions.checkNotNull(vcfAllele);
    if (m_chromosomeHgvsNameList == null) {
      m_chromosomeHgvsNameList = HGVS_NAME_SPLITTER.splitToList(m_chromosomeHgvsName);
    }
    if (vcfAllele.equals(".")) {
      // if phased, we should use "[?]"
      return "g." + m_position + "?";
    }
    if (vcfAllele.equals(m_ref)) {
      return "g." + m_position + "=";
    }
    for (int x = 0; x < m_alts.size(); x += 1) {
      if (vcfAllele.equals(m_alts.get(x))) {
        return m_chromosomeHgvsNameList.get(x);
      }
    }
    return "g." + m_position + m_ref + ">" + vcfAllele;
  }


  @Override
  public boolean equals(Object o) {
    if (this == o) {
      return true;
    }
    if (!(o instanceof VariantLocus that)) {
      return false;
    }
    return m_cpicPosition == that.getCpicPosition() &&
        m_position == that.getPosition() &&
        Objects.equals(m_chromosomeHgvsName, that.getChromosomeHgvsName()) &&
        Objects.equals(m_rsid, that.getRsid());
  }

  @Override
  public int hashCode() {
    return Objects.hash(m_cpicPosition, m_position, m_chromosomeHgvsName, m_rsid);
  }


  @Override
  public int compareTo(VariantLocus o) {

    int rez = ChromosomeNameComparator.getComparator().compare(m_chromosome, o.getChromosome());
    if (rez != 0) {
      return rez;
    }
    // WARNING: position can change (in DataManager), so any sorted collections in that code path must be careful
    rez = Long.compare(m_position, o.getPosition());
    if (rez != 0) {
      return rez;
    }
    rez = Long.compare(m_cpicPosition, o.getCpicPosition());
    if (rez != 0) {
      return rez;
    }
    return m_chromosomeHgvsName.compareTo(o.getChromosomeHgvsName());
  }

  @Override
  public String toString() {
    return String.format(
        "%s:%d%s",
        m_chromosome,
        m_cpicPosition,
        StringUtils.isBlank(m_rsid) ? "" : String.format(" (%s)", m_rsid));
  }
}

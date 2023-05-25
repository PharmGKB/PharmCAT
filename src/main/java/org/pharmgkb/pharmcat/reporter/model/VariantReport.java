package org.pharmgkb.pharmcat.reporter.model;

import java.lang.invoke.MethodHandles;
import java.util.Arrays;
import java.util.Collection;
import java.util.Set;
import java.util.TreeSet;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;
import org.apache.commons.lang3.ObjectUtils;
import org.apache.commons.lang3.StringUtils;
import org.pharmgkb.pharmcat.definition.model.VariantLocus;
import org.pharmgkb.pharmcat.haplotype.model.Variant;
import org.pharmgkb.pharmcat.util.HaplotypeNameComparator;
import org.pharmgkb.pharmcat.util.VariantUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


/**
 * A class to gather and standardize variant information needed for the final report.
 * <p>
 * Variant information can come from a few different sources so this class help give a unified interface for the report
 * to depend on.
 *
 * @author Ryan Whaley
 */
public class VariantReport implements Comparable<VariantReport> {
  private static final Logger sf_logger = LoggerFactory.getLogger(MethodHandles.lookup().lookupClass());
  
  @Expose
  @SerializedName("gene")
  private String m_gene;
  @Expose
  @SerializedName("chromosome")
  private String m_chr;
  @Expose
  @SerializedName("position")
  private long m_position;
  @Expose
  @SerializedName("dbSnpId")
  private String m_dbSnpId;
  @Expose
  @SerializedName("call")
  private String m_call;
  @Expose
  @SerializedName("alleles")
  private final Set<String> m_alleles = new TreeSet<>(HaplotypeNameComparator.getComparator());
  @Expose
  @SerializedName("phased")
  private boolean m_phased = false;
  @Expose
  @SerializedName("wildtypeAllele")
  private String m_wildTypeAllele;
  @Expose
  @SerializedName("hasUndocumentedVariations")
  private boolean m_hasUndocumentedVariations;
  @Expose
  @SerializedName("warnings")
  private Set<String> m_warnings = new TreeSet<>();

  public VariantReport(String gene, Variant variant) {
    setGene(gene);
    setPosition(variant.getPosition());
    setCall(variant.getVcfCall());
    setDbSnpId(variant.getRsid());
    setPhased(variant.isPhased());
  }

  public VariantReport(String gene, VariantLocus locus) {
    setGene(gene);
    setPosition(locus.getPosition());
    setDbSnpId(locus.getRsid());
  }

  public String getGene() {
    return m_gene;
  }

  public void setGene(String gene) {
    m_gene = gene;
  }

  public String getChr() {
    return m_chr;
  }

  public void setChr(String chr) {
    m_chr = chr;
  }

  public long getPosition() {
    return m_position;
  }

  public void setPosition(long position) {
    m_position = position;
  }

  public String getCall() {
    return m_call;
  }

  public void setCall(String call) {
    if (VariantUtils.isValidCall(call)) {
      m_call = call;
    }
    else {
      sf_logger.debug("Bad call value for {}: {}", this, call);
    }
  }

  public boolean isHetCall() {
    return VariantUtils.isHetCall(m_call);
  }

  public Set<String> getAlleles() {
    return m_alleles;
  }

  public void setAlleles(Collection<String> alleles) {
    m_alleles.addAll(alleles);
  }

  public boolean isPhased() {
    return m_phased;
  }

  public void setPhased(boolean phased) {
    m_phased = phased;
  }

  public String getDbSnpId() {
    return m_dbSnpId;
  }

  public void setDbSnpId(String dbSnpId) {
    m_dbSnpId = dbSnpId;
  }

  public String getWildTypeAllele() {
    return m_wildTypeAllele;
  }

  public void setWildTypeAllele(String wildTypeAllele) {
    m_wildTypeAllele = wildTypeAllele;
  }

  public boolean isMissing() {
    return StringUtils.isBlank(m_call);
  }
  
  public boolean isHasUndocumentedVariations() {
    return m_hasUndocumentedVariations;
  }
  
  public void setHasUndocumentedVariations(boolean hasUndocumentedVariations) {
    m_hasUndocumentedVariations = hasUndocumentedVariations;
  }

  public Set<String> getWarnings() {
    return m_warnings;
  }

  public void setWarnings(Collection<String> warnings) {
    m_warnings.clear();
    m_warnings.addAll(warnings);
  }

  public void addWarning(String warning) {
    m_warnings.add(warning);
  }

  public boolean isNonWildType() {
    return !(isMissing() || m_wildTypeAllele == null)
        && Arrays.stream(getCall().split("[|/]")).anyMatch(c -> !c.equals(getWildTypeAllele()));
  }


  @Override
  public String toString() {
    if (m_dbSnpId != null) {
      return m_dbSnpId;
    }
    return m_gene + ":" + m_position;
  }

  public String toChrPosition() {
    if (m_chr != null) {
      return m_chr + ":" + m_position;
    } else {
      return String.valueOf(m_position);
    }
  }

  @Override
  public int compareTo(VariantReport o) {
    // missing values get sorted last
    int rez = ObjectUtils.compare(isMissing(), o.isMissing());
    if (rez != 0) {
      return rez;
    }
    
    // then order by chromosome
    rez = ObjectUtils.compare(getChr(), o.getChr());
    if (rez != 0) {
      return rez;
    }

    // then by position on the chromosome
    return ObjectUtils.compare(getPosition(), o.getPosition());
  }
}

package org.pharmgkb.pharmcat.reporter.model;

import java.util.Collection;
import java.util.Set;
import java.util.TreeSet;
import javax.annotation.Nonnull;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;
import org.apache.commons.lang3.StringUtils;
import org.pharmgkb.common.comparator.HaplotypeNameComparator;
import org.pharmgkb.pharmcat.definition.VariantAlleleMap;
import org.pharmgkb.pharmcat.definition.model.VariantLocus;
import org.pharmgkb.pharmcat.haplotype.model.Variant;


/**
 * A class to gather and standardize variant information needed for the final report
 *
 * Variant information can come from a few different sources so this class help give a unified interface for the report
 * to depend on.
 *
 * @author Ryan Whaley
 */
public class VariantReport {

  @Expose
  @SerializedName("gene")
  private String m_gene;
  @Expose
  @SerializedName("position")
  private int m_position;
  @Expose
  @SerializedName("dbSnpId")
  private String m_dbSnpId;
  @Expose
  @SerializedName("call")
  private String m_call;
  @Expose
  @SerializedName("alleles")
  private Set<String> m_alleles = new TreeSet<>(HaplotypeNameComparator.getComparator());
  @Expose
  @SerializedName("phased")
  private boolean m_phased = false;

  public VariantReport(String gene, Variant variant) {
    m_gene = gene;
    m_position = variant.getPosition();
    m_call = variant.getVcfCall();
    m_dbSnpId = variant.getRsid();
    m_phased = variant.isPhased();
  }

  public VariantReport(String gene, VariantLocus locus) {
    m_gene = gene;
    m_position = locus.getPosition();
    m_dbSnpId = locus.getRsid();
  }

  public VariantReport findAlleles(@Nonnull VariantAlleleMap variantAlleleMap) {
    Collection<String> alleles = variantAlleleMap.getAlleles(m_position);

    if (alleles != null) {
      m_alleles.addAll(alleles);
    }

    return this;
  }

  public String getGene() {
    return m_gene;
  }

  public void setGene(String gene) {
    m_gene = gene;
  }

  public int getPosition() {
    return m_position;
  }

  public void setPosition(int position) {
    m_position = position;
  }

  public String getCall() {
    return m_call;
  }

  public void setCall(String call) {
    m_call = call;
  }

  public Set<String> getAlleles() {
    return m_alleles;
  }

  public void setAlleles(Set<String> alleles) {
    m_alleles = alleles;
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

  public boolean isMissing() {
    return StringUtils.isBlank(m_call);
  }
}

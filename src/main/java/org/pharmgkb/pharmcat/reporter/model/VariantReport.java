package org.pharmgkb.pharmcat.reporter.model;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;
import org.apache.commons.lang3.StringUtils;
import org.pharmgkb.common.comparator.HaplotypeNameComparator;
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
  @Expose
  @SerializedName("wildtypeAllele")
  private String m_wildtypeAllele;
  @Expose
  @SerializedName("messages")
  private List<MessageVariant> m_messages = new ArrayList<>();

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

  public String getWildtypeAllele() {
    return m_wildtypeAllele;
  }

  public void setWildtypeAllele(String wildtypeAllele) {
    m_wildtypeAllele = wildtypeAllele;
  }

  public List<MessageVariant> getMessages() {
    return m_messages;
  }

  public void addMessage(MessageVariant message) {
    m_messages.add(message);
  }

  public boolean isMissing() {
    return StringUtils.isBlank(m_call);
  }

  public boolean isNonwildtype() {
    return !(isMissing() || m_wildtypeAllele == null)
        && Arrays.stream(getCall().split("[|/]")).anyMatch(c -> !c.equals(getWildtypeAllele()));
  }

  @Override
  public String toString() {
    if (m_dbSnpId != null) {
      return m_dbSnpId;
    }
    return m_gene + ":" + m_position;
  }
}

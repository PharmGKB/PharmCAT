package org.pharmgkb.pharmcat.reporter.model;

import java.util.ArrayList;
import java.util.List;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;
import org.apache.commons.lang3.StringUtils;
import org.jspecify.annotations.Nullable;
import org.pharmgkb.common.util.ComparisonChain;


/**
 * @author Ryan Whaley
 */
public class MatchLogic implements Comparable<MatchLogic> {

  @Expose
  @SerializedName("gene")
  private String m_gene;
  @Expose
  @SerializedName("hapsCalled")
  private final List<String> m_hapsCalled = new ArrayList<>();
  @Expose
  @SerializedName("hapsMissing")
  private final List<String> m_hapsMissing = new ArrayList<>();
  @Expose
  @SerializedName("variantsMissing")
  private final List<String> m_variantsMissing = new ArrayList<>();
  @Expose
  @SerializedName("variant")
  private List<String> m_variants = new ArrayList<>();
  @Expose
  @SerializedName("dips")
  private final List<String> m_dips = new ArrayList<>();
  @Expose
  @SerializedName("drugs")
  private final List<String> m_drugs = new ArrayList<>();


  public @Nullable String getGene() {
    return m_gene;
  }

  public void setGene(@Nullable String gene) {
    m_gene = StringUtils.stripToNull(gene);
  }

  public List<String> getHapsCalled() {
    return m_hapsCalled;
  }

  public void setHapsCalled(List<String> haps) {
    m_hapsCalled.addAll(haps);
  }

  public List<String> getDips() {
    return m_dips;
  }

  public void setDips(List<String> dips) {
    m_dips.addAll(dips);
  }

  public List<String> getDrugs() {
    return m_drugs;
  }

  public void setDrugs(List<String> drugs) {
    m_drugs.addAll(drugs);
  }

  public List<String> getHapsMissing() {
    return m_hapsMissing;
  }

  public void setHapsMissing(List<String> hapsMissing) {
    m_hapsMissing.addAll(hapsMissing);
  }

  public List<String> getVariantsMissing() {
    return m_variantsMissing;
  }

  public void setVariantsMissing(List<String> variantsMissing) {
    m_variantsMissing.addAll(variantsMissing);
  }

  public List<String> getVariants() {
    return m_variants;
  }

  public void setVariants(List<String> variants) {
    m_variants = variants;
  }


  @Override
  public int compareTo(MatchLogic o) {
    if (o == this) {
      return 0;
    }
    return new ComparisonChain()
        .compare(m_gene, o.getGene())
        .compare(m_hapsCalled, o.getHapsCalled())
        .compare(m_hapsMissing, o.getHapsMissing())
        .compare(m_dips, o.getDips())
        .compare(m_drugs, o.getDrugs())
        .compare(m_variants, o.getVariants())
        .compare(m_variantsMissing, o.getVariantsMissing())
        .result();
  }
}

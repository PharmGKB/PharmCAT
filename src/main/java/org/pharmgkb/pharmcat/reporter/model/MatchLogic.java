package org.pharmgkb.pharmcat.reporter.model;

import java.util.ArrayList;
import java.util.List;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;
import org.apache.commons.lang3.StringUtils;
import org.checkerframework.checker.nullness.qual.Nullable;


/**
 * @author Ryan Whaley
 */
public class MatchLogic {

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
  private String m_variant;
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

  public @Nullable String getVariant() {
    return m_variant;
  }

  public void setVariant(@Nullable String variant) {
    m_variant = StringUtils.stripToNull(variant);
  }
}

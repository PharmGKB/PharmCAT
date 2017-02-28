package org.pharmgkb.pharmcat.reporter.model;

import java.util.ArrayList;
import java.util.List;
import javax.annotation.Nonnull;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;


/**
 * @author Ryan Whaley
 */
public class MatchLogic {
  
  @Expose
  @SerializedName("gene")
  private String m_gene;
  @Expose
  @SerializedName("hapsCalled")
  private List<String> m_hapsCalled = new ArrayList<>();
  @Expose
  @SerializedName("hapsMissing")
  private List<String> m_hapsMissing = new ArrayList<>();
  @Expose
  @SerializedName("variantsMissing")
  private List<String> m_variantsMissing = new ArrayList<>();
  @Expose
  @SerializedName("dips")
  private List<String> m_dips = new ArrayList<>();
  @Expose
  @SerializedName("hapsAvailable")
  private List<String> m_hapsAvailable = new ArrayList<>();
  @Expose
  @SerializedName("drugs")
  private List<String> m_drugs = new ArrayList<>();


  @Nonnull
  public String getGene() {
    return m_gene;
  }

  public void setGene(String gene) {
    m_gene = gene;
  }

  @Nonnull
  public List<String> getHapsCalled() {
    return m_hapsCalled;
  }

  public void setHapsCalled(List<String> haps) {
    m_hapsCalled.addAll(haps);
  }

  @Nonnull
  public List<String> getHapsAvailable() {
    return m_hapsAvailable;
  }

  public void setHapsAvailable(List<String> haps) {
    m_hapsAvailable.addAll(haps);
  }

  @Nonnull
  public List<String> getDips() {
    return m_dips;
  }

  public void setDips(List<String> dips) {
    m_dips.addAll(dips);
  }

  @Nonnull
  public List<String> getDrugs() {
    return m_drugs;
  }

  public void setDrugs(List<String> drugs) {
    m_drugs = drugs;
  }

  public List<String> getHapsMissing() {
    return m_hapsMissing;
  }

  public void setHapsMissing(List<String> hapsMissing) {
    m_hapsMissing = hapsMissing;
  }

  public List<String> getVariantsMissing() {
    return m_variantsMissing;
  }

  public void setVariantsMissing(List<String> variantsMissing) {
    m_variantsMissing = variantsMissing;
  }
}

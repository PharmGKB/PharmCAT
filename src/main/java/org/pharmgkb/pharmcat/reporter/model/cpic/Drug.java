package org.pharmgkb.pharmcat.reporter.model.cpic;

import java.util.List;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;


/**
 * A CPIC Drug with related information sourced from the CPIC database.
 *
 * This includes:
 * <ul>
 *   <li>recommendations</li>
 *   <li>citations</li>
 *   <li>guideline info</li>
 *   <li>related genes</li>
 * </ul>
 */
public class Drug {
  @Expose
  @SerializedName("drugid")
  private String m_drugId;
  @Expose
  @SerializedName("drugname")
  private String m_drugName;
  @Expose
  @SerializedName("guidelinename")
  private String m_guidelineName;
  @Expose
  @SerializedName("url")
  private String m_url;
  @Expose
  @SerializedName("guidelinepharmgkbids")
  private List<String> m_guidelinePharmgkbIds;
  @Expose
  @SerializedName("citations")
  private List<Publication> m_publications;
  @Expose
  @SerializedName("recommendations")
  private List<Recommendation> m_recommendations;
  @Expose
  @SerializedName("genes")
  private List<String> m_genes;

  public String toString() {
    return m_drugName;
  }

  public String getDrugId() {
    return m_drugId;
  }

  public void setDrugId(String drugId) {
    m_drugId = drugId;
  }

  public String getDrugName() {
    return m_drugName;
  }

  public void setDrugName(String drugName) {
    m_drugName = drugName;
  }

  public String getGuidelineName() {
    return m_guidelineName;
  }

  public void setGuidelineName(String guidelineName) {
    m_guidelineName = guidelineName;
  }

  public String getUrl() {
    return m_url;
  }

  public void setUrl(String url) {
    m_url = url;
  }

  public List<String> getGuidelinePharmgkbIds() {
    return m_guidelinePharmgkbIds;
  }

  public void setGuidelinePharmgkbIds(List<String> guidelinePharmgkbIds) {
    m_guidelinePharmgkbIds = guidelinePharmgkbIds;
  }

  public List<Publication> getPublications() {
    return m_publications;
  }

  public void setPublications(List<Publication> publications) {
    m_publications = publications;
  }

  public List<Recommendation> getRecommendations() {
    return m_recommendations;
  }

  public void setRecommendations(List<Recommendation> recommendations) {
    m_recommendations = recommendations;
  }

  public List<String> getGenes() {
    return m_genes;
  }

  public void setGenes(List<String> genes) {
    m_genes = genes;
  }
}

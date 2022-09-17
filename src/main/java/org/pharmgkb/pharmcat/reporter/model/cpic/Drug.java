package org.pharmgkb.pharmcat.reporter.model.cpic;

import java.util.List;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;
import org.apache.commons.lang3.ObjectUtils;
import org.pharmgkb.pharmcat.reporter.model.DataSource;


/**
 * A drug record from an outside data source with related information.
 *
 * This includes:
 * <ul>
 *   <li>recommendations</li>
 *   <li>citations</li>
 *   <li>guideline info</li>
 *   <li>related genes</li>
 * </ul>
 */
public class Drug implements Comparable<Drug> {
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
  private List<Publication> m_citations;
  @Expose
  @SerializedName("recommendations")
  private List<Recommendation> m_recommendations;
  @Expose
  @SerializedName("genes")
  private List<String> m_genes;
  @Expose
  @SerializedName("notesonusage")
  private String m_notesOnUsage;
  @Expose
  @SerializedName("cpicVersion")
  private String m_cpicVersion;
  @Expose
  @SerializedName("source")
  private DataSource m_source;

  public String toString() {
    return String.format("%s [%s]", m_drugName, m_source);
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

  public List<Publication> getCitations() {
    return m_citations;
  }

  public void setCitations(List<Publication> citations) {
    m_citations = citations;
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

  public String getNotesOnUsage() {
    return m_notesOnUsage;
  }

  public void setNotesOnUsage(String notesOnUsage) {
    m_notesOnUsage = notesOnUsage;
  }

  public String getCpicVersion() {
    return m_cpicVersion;
  }

  public void setCpicVersion(String cpicVersion) {
    m_cpicVersion = cpicVersion;
  }

  public DataSource getSource() {
    return m_source;
  }

  public void setSource(DataSource source) {
    m_source = source;

    if (m_recommendations != null) {
      m_recommendations.forEach(r -> r.setSource(source));
    }
  }

  @Override
  public int compareTo(Drug other) {
    int rez = m_drugName.compareToIgnoreCase(other.getDrugName());
    if (rez != 0) {
      return rez;
    }
    rez = ObjectUtils.compare(m_source, other.getSource());
    if (rez != 0) {
      return rez;
    }
    return m_guidelineName.compareToIgnoreCase(other.getGuidelineName());
  }
}

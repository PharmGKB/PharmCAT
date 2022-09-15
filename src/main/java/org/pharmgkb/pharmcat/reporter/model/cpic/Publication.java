package org.pharmgkb.pharmcat.reporter.model.cpic;

import java.util.List;
import java.util.Map;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;


/**
 * A Publication object sourced from the CPIC DB
 */
public class Publication {
  @Expose
  @SerializedName("pmid")
  private String m_pmid;
  @Expose
  @SerializedName("title")
  private String m_title;
  @Expose
  @SerializedName("authors")
  private List<String> m_authors;
  @Expose
  @SerializedName("journal")
  private String m_journal;
  @Expose
  @SerializedName("year")
  private Long m_year;
  // this only comes in from PharmGKB data
  @Expose(serialize = false)
  @SerializedName("crossReferences")
  private List<Map<String, Object>> m_crossReferences;


  public String getPmid() {
    return m_pmid;
  }

  public String getTitle() {
    return m_title;
  }

  public List<String> getAuthors() {
    return m_authors;
  }

  public String getJournal() {
    return m_journal;
  }

  public Long getYear() {
    return m_year;
  }


  /**
   * This massages the original data for use in PharmCAT.
   */
  public void normalize() {
    if (m_pmid == null && m_crossReferences != null) {
      for (Map<String, Object> xref : m_crossReferences) {
        if (xref.get("resource").equals("PubMed")) {
          m_pmid = (String)xref.get("resourceId");
          break;
        }
      }
    }
  }
}

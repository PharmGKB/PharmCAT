package org.pharmgkb.pharmcat.reporter.model;

import java.util.List;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;


/**
 * Literature object used for citations
 *
 * @author Ryan Whaley
 */
public class Literature {

  @SerializedName("title")
  @Expose
  private String m_title;
  @SerializedName("authors")
  @Expose
  private List<String> m_authors;
  @SerializedName("journal")
  @Expose
  private String m_journal;
  @SerializedName("year")
  @Expose
  private Integer m_year;
  @SerializedName("resourceId")
  @Expose
  private String m_pmid;

  public String getTitle() {
    return m_title;
  }

  public List<String> getAuthors() {
    return m_authors;
  }

  public String getJournal() {
    return m_journal;
  }

  public Integer getYear() {
    return m_year;
  }

  public String getPmid() {
    return m_pmid;
  }
}

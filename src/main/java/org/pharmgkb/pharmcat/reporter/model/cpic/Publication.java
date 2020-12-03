package org.pharmgkb.pharmcat.reporter.model.cpic;

import java.util.List;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;


/**
 * A Publication object sourced from the CPIC DB
 */
public class Publication {
  @Expose
  @SerializedName("journal")
  private String m_journal;
  @Expose
  @SerializedName("title")
  private String m_title;
  @Expose
  @SerializedName("authors")
  private List<String> m_auhtors;
  @Expose
  @SerializedName("pmid")
  private String m_pmid;
  @Expose
  @SerializedName("pmcid")
  private String m_pmcid;
  @Expose
  @SerializedName("doi")
  private String m_doi;
  @Expose
  @SerializedName("month")
  private Long m_month;
  @Expose
  @SerializedName("year")
  private Long m_year;

  public String getJournal() {
    return m_journal;
  }

  public void setJournal(String journal) {
    m_journal = journal;
  }

  public String getTitle() {
    return m_title;
  }

  public void setTitle(String title) {
    m_title = title;
  }

  public List<String> getAuhtors() {
    return m_auhtors;
  }

  public void setAuhtors(List<String> auhtors) {
    m_auhtors = auhtors;
  }

  public String getPmid() {
    return m_pmid;
  }

  public void setPmid(String pmid) {
    m_pmid = pmid;
  }

  public String getPmcid() {
    return m_pmcid;
  }

  public void setPmcid(String pmcid) {
    m_pmcid = pmcid;
  }

  public String getDoi() {
    return m_doi;
  }

  public void setDoi(String doi) {
    m_doi = doi;
  }

  public Long getMonth() {
    return m_month;
  }

  public void setMonth(Long month) {
    m_month = month;
  }

  public Long getYear() {
    return m_year;
  }

  public void setYear(Long year) {
    m_year = year;
  }
}

package org.pharmgkb.pharmcat.reporter.model.result;

import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;


/**
 * @author Ryan Whaley
 * @deprecated
 */
public class Markdown {
  
  @SerializedName("id")
  @Expose
  private long m_id = -1;
  @SerializedName("html")
  @Expose
  private String m_html;

  public String getHtml() {
    return m_html;
  }

  public void setHtml(String html) {
    m_html = html;
  }

  public long getId() {
    return m_id;
  }

  public void setId(long id) {
    m_id = id;
  }
}

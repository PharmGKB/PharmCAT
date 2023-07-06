package org.pharmgkb.pharmcat.reporter.model.pgkb;

import java.util.regex.Matcher;
import java.util.regex.Pattern;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;
import org.apache.commons.lang3.StringUtils;


/**
 * PharmGKB Markdown Model
 * @author Ryan Whaley
 */
public class Markdown {
  private static final Pattern PTAG_PATTERN = Pattern.compile("</?p>");

  @SerializedName("id")
  @Expose
  private long m_id = -1;
  @SerializedName("html")
  @Expose
  private String m_html;


  public long getId() {
    return m_id;
  }

  public String getHtml() {
    return m_html;
  }

  /**
   * Strips the opening and closing <code>p</code> tags from an HTML String. Any inner markup will remain.
   * @return a stripped HTML String
   */
  public String getHtmlStripped() {
    Matcher m = PTAG_PATTERN.matcher(StringUtils.strip(m_html));
    return m.replaceAll("");
  }
}

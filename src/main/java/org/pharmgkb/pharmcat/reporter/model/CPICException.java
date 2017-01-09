package org.pharmgkb.pharmcat.reporter.model;

import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;


/**
 * @author Lester Carter
 */
public class CPICException {

  @Expose
  @SerializedName("rule_name")
  private String m_name;
  @Expose
  @SerializedName("version")
  private String m_version;
  @Expose
  @SerializedName("matches")
  private CPICExceptionMatch m_matches;
  @Expose
  @SerializedName("exception_type")
  private String m_exceptionType;
  @Expose
  @SerializedName("message")
  private String m_message;

  public String getName() {
    return m_name;
  }

  public void setName(String rule_name) {
    m_name = rule_name;
  }

  public String getVersion() {
    return m_version;
  }

  public void setVersion(String version) {
    m_version = version;
  }

  public CPICExceptionMatch getMatches() {
    return m_matches;
  }

  public void setMatches(CPICExceptionMatch matches) {
    m_matches = matches;
  }

  public String getExceptionType() {
    return m_exceptionType;
  }

  public void setExceptionType(String exceptionType) {
    m_exceptionType = exceptionType;
  }

  public String getMessage() {
    return m_message;
  }

  public void setMessage(String message) {
    m_message = message;
  }
}

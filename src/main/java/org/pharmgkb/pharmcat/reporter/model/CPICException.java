package org.pharmgkb.pharmcat.reporter.model;

import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;


/**
 * @author Lester Carter
 */
public class CPICException {

  @SerializedName("rule_name")
  @Expose
  private String m_name;
  String version;
  String gene;
  String matches;
  String exception_type;
  String message;

  public String getName() {
    return m_name;
  }

  public void setName(String rule_name) {
    m_name = rule_name;
  }

  public String getVersion() {
    return version;
  }

  public void setVersion(String version) {
    this.version = version;
  }

  public String getGene() {
    return gene;
  }

  public void setGene(String gene) {
    this.gene = gene;
  }

  public String getMatches() {
    return matches;
  }

  public void setMatches(String matches) {
    this.matches = matches;
  }

  public String getException_type() {
    return exception_type;
  }

  public void setException_type(String exception_type) {
    this.exception_type = exception_type;
  }

  public String getMessage() {
    return message;
  }

  public void setMessage(String message) {
    this.message = message;
  }




}

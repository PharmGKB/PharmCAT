package org.cpic.reporter.model;

/**
 * Created by lester on 4/13/16.
 */
public class CPICException {

  String rule_name;
  String version;
  String gene;
  String matches;
  String exception_type;
  String message;

  public String getRule_name() {
    return rule_name;
  }

  public void setRule_name(String rule_name) {
    this.rule_name = rule_name;
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

package org.pharmgkb.pharmcat.reporter.model;

import java.util.ArrayList;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;


/**
 * @author Lester Carter
 */
public class CPICExceptionList {

  @SerializedName("rules")
  @Expose
  private ArrayList<CPICException> m_rules;

  public ArrayList<CPICException> getRules() {
    return m_rules;
  }

  public void setRules(ArrayList<CPICException> rules) {
    m_rules = rules;
  }
}

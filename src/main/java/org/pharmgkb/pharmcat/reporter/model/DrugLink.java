package org.pharmgkb.pharmcat.reporter.model;

import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;


/**
 * Model object to hold all properties needed to display a link to a drug
 *
 * @author Ryan Whaley
 */
public class DrugLink {

  @Expose
  @SerializedName("name")
  private String m_name;

  @Expose
  @SerializedName("guidelineId")
  private String m_guidelineId;

  public DrugLink(String name, String guidelineId) {
    m_name = name;
    m_guidelineId = guidelineId;
  }

  public String getName() {
    return m_name;
  }

  public String getGuidelineId() {
    return m_guidelineId;
  }
}

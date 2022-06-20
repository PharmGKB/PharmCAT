package org.pharmgkb.pharmcat.reporter.model;

import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;
import org.apache.commons.lang3.builder.EqualsBuilder;


/**
 * Model object to hold all properties needed to display a link to a drug
 *
 * @author Ryan Whaley
 */
public class DrugLink implements Comparable<DrugLink> {

  @Expose
  @SerializedName("name")
  private String m_name;

  @Expose
  @SerializedName("id")
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

  public boolean equals(Object o) {
    if (!(o instanceof DrugLink)) {
      return false;
    }
    if (o == this) {
      return true;
    }

    DrugLink d = (DrugLink)o;
    return new EqualsBuilder()
        .append(m_guidelineId, d.getGuidelineId())
        .append(m_name, d.getName())
        .isEquals();
  }

  public int compareTo(DrugLink o) {
    if (o == null) {
      return -1;
    }
    int rez = m_name.compareTo(o.getName());

    if (rez != 0) {
      return rez;
    } else {
      return m_guidelineId.compareTo(o.getGuidelineId());
    }
  }
}

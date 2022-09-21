package org.pharmgkb.pharmcat.reporter.model.result;

import java.util.Objects;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;
import org.apache.commons.lang3.builder.CompareToBuilder;


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


  @Override
  public int compareTo(DrugLink o) {
    return new CompareToBuilder()
        .append(m_name, o.getName())
        .append(m_guidelineId, o.getGuidelineId())
        .toComparison();
  }

  @Override
  public boolean equals(Object o) {
    if (!(o instanceof DrugLink d)) {
      return false;
    }
    if (o == this) {
      return true;
    }

    return Objects.equals(m_name, d.getName()) &&
        Objects.equals(m_guidelineId, d.getGuidelineId());
  }

  @Override
  public int hashCode() {
    return Objects.hash(m_name, m_guidelineId);
  }
}

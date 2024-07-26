package org.pharmgkb.pharmcat.reporter.model.pgkb;

import java.util.List;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;


public class PrescribingGuidanceDataset {
  @Expose
  @SerializedName("version")
  private String m_version;

  @Expose
  @SerializedName("guidelines")
  private List<GuidelinePackage> m_guidelinePackages;


  public String getVersion() {
    return m_version;
  }

  public List<GuidelinePackage> getGuidelinePackages() {
    return m_guidelinePackages;
  }

  public String toString() {
    return "Prescribing Guidance " + getVersion();
  }
}

package org.pharmgkb.pharmcat.reporter.model;

import java.util.ArrayList;
import java.util.List;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;


/**
 * This class packs together the main guideline and the groups of diplotype-specific annotations
 *
 * @author Ryan Whaley
 */
public class GuidelinePackage {

  @SerializedName("guideline")
  @Expose
  private DosingGuideline guideline;
  @SerializedName("annotationGroups")
  @Expose
  private List<Group> groups = new ArrayList<>();

  public DosingGuideline getGuideline() {
    return guideline;
  }

  public void setGuideline(DosingGuideline guideline) {
    this.guideline = guideline;
  }

  public List<Group> getGroups() {
    return groups;
  }

  public void setGroups(List<Group> groups) {
    this.groups = groups;
  }
}

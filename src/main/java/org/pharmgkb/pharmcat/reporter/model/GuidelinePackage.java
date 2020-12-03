package org.pharmgkb.pharmcat.reporter.model;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;


/**
 * This class packs together the main guideline and the groups of diplotype-specific annotations
 *
 * @author Ryan Whaley
 * @deprecated
 */
public class GuidelinePackage {

  @SerializedName("guideline")
  @Expose
  private DosingGuideline guideline;
  @SerializedName("annotationGroups")
  @Expose
  private List<Group> groups = new ArrayList<>();
  @SerializedName("genePhenotypes")
  @Expose
  private Map<String,Map<String,String>> phenotypeMap = new HashMap<>();
  @SerializedName("citations")
  @Expose
  private List<Literature> citations = new ArrayList<>();

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

  public List<Literature> getCitations() {
    return citations;
  }

  /**
   * This is a map of: "Gene Symbol" -&gt; "Allele Name" -&gt; "Phenotype". This can be used in conjuction with the
   * "genePhenotypes" property in the {@link Group} object.
   * @return a Mapping of gene+allele -&gt; phenotype
   */
  public Map<String, Map<String,String>> getPhenotypeMap() {
    return phenotypeMap;
  }

  public void setPhenotypeMap(Map<String, Map<String,String>> phenotypeMap) {
    this.phenotypeMap = phenotypeMap;
  }

  @Override
  public String toString() {
    if (guideline != null) {
      return guideline.getName();
    }
    else {
      return super.toString();
    }
  }
}

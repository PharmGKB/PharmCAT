package org.pharmgkb.pharmcat.reporter.model.pgkb;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;
import org.pharmgkb.pharmcat.reporter.model.cpic.Publication;
import org.pharmgkb.pharmcat.reporter.model.result.Genotype;


/**
 * This class packs together the main Guideline Annotation and the Groups of diplotype-specific annotations
 *
 * @author Ryan Whaley
 */
public class GuidelinePackage implements Comparable<GuidelinePackage> {
  @Expose
  @SerializedName("guideline")
  private DosingGuideline guideline;
  @Expose
  @SerializedName("annotationGroups")
  private List<Group> groups = new ArrayList<>();
  @Expose
  @SerializedName("citations")
  private List<Publication> citations = new ArrayList<>();
  @Expose
  @SerializedName("version")
  private String m_version;


  /**
   * Private constructor for GSON;
   */
  private GuidelinePackage() {
  }


  public DosingGuideline getGuideline() {
    return guideline;
  }


  public List<Group> getGroups() {
    return groups;
  }


  public List<Publication> getCitations() {
    return citations;
  }

  
  public Set<String> getGenes() {
    return guideline.getGuidelineGenes().stream()
        .map(g -> g.getGene().getSymbol())
        .collect(Collectors.toSet());
  }

  public Set<String> getDrugs() {
    return guideline.getRelatedChemicals().stream()
        .map(AccessionObject::getName)
        .collect(Collectors.toSet());
  }


  public String getVersion() {
    return m_version;
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

  @Override
  public int compareTo(GuidelinePackage o) {
    return toString().compareToIgnoreCase(o.toString());
  }


  public void applyFunctions(Genotype genotype) {
    getGuideline().applyFunctions(genotype);
  }
}

package org.pharmgkb.pharmcat.reporter.model.pgkb;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;
import org.pharmgkb.pharmcat.reporter.model.result.Diplotype;


/**
 * This class packs together the main Guideline Annotation and the Groups of diplotype-specific annotations
 *
 * @author Ryan Whaley
 */
public class GuidelinePackage implements Comparable<GuidelinePackage> {

  @SerializedName("guideline")
  @Expose
  private DosingGuideline guideline;
  @SerializedName("annotationGroups")
  @Expose
  private List<Group> groups = new ArrayList<>();
  @SerializedName("matchedGroups")
  @Expose
  private List<Group> matchedGroups = new ArrayList<>();
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

  
  public List<Group> getMatchedGroups() {
    return matchedGroups;
  }

  
  public List<Literature> getCitations() {
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

  
  public void matchGroups(Collection<Diplotype> diplotypes) {
    for (Diplotype diplotype : diplotypes) {
      if (!diplotype.isUnknownPhenotype() && diplotype.isUnknownAlleles()) {
        getGroups().stream()
            .filter(group -> diplotype.getPhenotypes().stream().anyMatch(p -> group.getName().equalsIgnoreCase(p)))
            .forEach(g -> {
              g.addMatchingDiplotype(diplotype);
              matchedGroups.add(g);
            });
      } else if (!diplotype.isUnknownAlleles()) {
        Set<String> functionKeys = getGuideline().getFunctionKeysForDiplotype(diplotype);
        for (String functionKey : functionKeys) {
          getGroups().stream()
              .filter(group -> group.matchesKey(functionKey))
              .forEach(g -> {
                g.addMatchingFunctionKey(functionKey);
                g.addMatchingDiplotype(diplotype);
                matchedGroups.add(g);
              });
        }
      }
    }
  }

  public boolean hasMatch() {
    return matchedGroups.size() > 0;
  }
}

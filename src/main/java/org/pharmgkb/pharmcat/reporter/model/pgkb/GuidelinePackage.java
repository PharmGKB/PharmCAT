package org.pharmgkb.pharmcat.reporter.model.pgkb;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;
import java.util.stream.Collectors;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;
import org.pharmgkb.pharmcat.reporter.model.cpic.Publication;
import org.pharmgkb.pharmcat.reporter.model.result.Diplotype;
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

  private final Set<Group> matchedGroups = new TreeSet<>();


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

  
  public Set<Group> getMatchedGroups() {
    return matchedGroups;
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

  
  public void match(Genotype genotype) {
    for (Diplotype diplotype : genotype.getDiplotypes()) {
      if (diplotype.isPhenotypeOnly() || diplotype.isAllelePresenceType()) {
        groups.stream()
            .filter(group -> diplotype.getPhenotypes().stream().anyMatch(p -> group.getName().equalsIgnoreCase(p)))
            .forEach(g -> {
              g.addMatchingDiplotype(diplotype);
              g.addMatchingGenotype(genotype);
              matchedGroups.add(g);
            });
      } else if (!diplotype.isUnknownAlleles()) {
        Set<String> functionKeys = guideline.getFunctionKeysForDiplotype(diplotype);
        for (String functionKey : functionKeys) {
          getGroups().stream()
              .filter(group -> group.matchesKey(functionKey))
              .forEach(g -> {
                g.addMatchingFunctionKey(functionKey);
                g.addMatchingDiplotype(diplotype);
                g.addMatchingGenotype(genotype);
                matchedGroups.add(g);
              });
        }
      }
    }
  }

  public void applyFunctions(Genotype genotype) {
    getGuideline().applyFunctions(genotype);
  }

  public boolean hasMatch() {
    return matchedGroups.size() > 0;
  }

  public Integer getVersion() {
    return guideline.getVersion();
  }
}

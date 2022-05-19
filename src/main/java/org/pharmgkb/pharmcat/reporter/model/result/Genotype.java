package org.pharmgkb.pharmcat.reporter.model.result;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.stream.Collectors;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;


public class Genotype {

  @Expose
  @SerializedName("diplotypes")
  private final SortedSet<Diplotype> diplotypes = new TreeSet<>();

  public Genotype(Collection<Diplotype> diplotypes) {
    diplotypes.forEach(this::addDiplotype);
  }

  public Genotype(Diplotype diplotype) {
    addDiplotype(diplotype);
  }

  public SortedSet<Diplotype> getDiplotypes() {
    return diplotypes;
  }

  public void addDiplotype(Diplotype diplotype) {
    if (getDiplotypes().stream().anyMatch(d -> d.getGene().equals(diplotype.getGene()))) {
      throw new RuntimeException("Genotype already has diplotype for " + diplotype.getGene());
    }
    this.diplotypes.add(diplotype);
  }

  public String toString() {
    if (this.diplotypes.size() == 0) {
      return "Empty Genotype";
    } else {
      return this.diplotypes.stream()
          .map(Diplotype::toString)
          .collect(Collectors.joining("; "));
    }
  }

  public List<Map<String,String>> toLookupKeys() {
    List<Map<String,String>> keys = new ArrayList<>();
    for (Diplotype diplotype : getDiplotypes()) {
      keys = addToLookupKeyList(diplotype, keys);
    }
    return keys;
  }

  private List<Map<String,String>> addToLookupKeyList(Diplotype diplotype, List<Map<String,String>> originalList) {
    List<Map<String,String>> newList = new ArrayList<>();
    for (String lookupKey : diplotype.getLookupKeys()) {
      if (originalList != null && originalList.size() > 0) {
        for (Map<String, String> originalKey : originalList) {
          Map<String,String> newLookupMap = new HashMap<>();
          newLookupMap.put(diplotype.getGene(), lookupKey);
          newLookupMap.putAll(originalKey);
          newList.add(newLookupMap);
        }
      } else {
        Map<String,String> newLookupMap = new HashMap<>();
        newLookupMap.put(diplotype.getGene(), lookupKey);
        newList.add(newLookupMap);
      }
    }
    return newList;
  }
}

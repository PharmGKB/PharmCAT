package org.pharmgkb.pharmcat.reporter.model.pgkb;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import javax.annotation.Nonnull;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;
import org.pharmgkb.pharmcat.reporter.RecommendationUtils;
import org.pharmgkb.pharmcat.reporter.model.result.Genotype;


public class RecommendationAnnotation implements Comparable<RecommendationAnnotation> {

  @Expose
  @SerializedName("id")
  private String id;
  @Expose
  @SerializedName("name")
  private String name;
  @Expose
  @SerializedName("population")
  private String population;
  @Expose
  @SerializedName("classification")
  private OntologyTerm classification;
  @Expose
  @SerializedName("relatedChemicals")
  private List<AccessionObject> relatedChemicals;
  @SerializedName("text")
  @Expose
  private Markdown text;
  @SerializedName("implications")
  @Expose
  private List<String> implications = new ArrayList<>();
  @SerializedName("lookupKey")
  @Expose
  private Map<String,Object> lookupKey = new HashMap<>();
  @SerializedName("dosingInformation")
  @Expose
  private boolean dosingInformation;
  @SerializedName("alternateDrugAvailable")
  @Expose
  private boolean alternateDrugAvailable;
  @SerializedName("otherPrescribingGuidance")
  @Expose
  private boolean otherPrescribingGuidance;


  public String getId() {
    return id;
  }

  public void setId(String id) {
    this.id = id;
  }


  public String getName() {
    return name;
  }

  public void setName(String name) {
    this.name = name;
  }


  public String getPopulation() {
    return population;
  }

  public void setPopulation(String population) {
    this.population = population;
  }


  public List<AccessionObject> getRelatedChemicals() {
    return relatedChemicals;
  }

  public void setRelatedChemicals(List<AccessionObject> relatedChemicals) {
    this.relatedChemicals = relatedChemicals;
  }


  public OntologyTerm getClassification() {
    return classification;
  }

  public void setClassification(OntologyTerm classification) {
    this.classification = classification;
  }


  public boolean isAlternateDrugAvailable() {
    return alternateDrugAvailable;
  }

  public void setAlternateDrugAvailable(boolean alternateDrugAvailable) {
    this.alternateDrugAvailable = alternateDrugAvailable;
  }


  public boolean isDosingInformation() {
    return dosingInformation;
  }

  public void setDosingInformation(boolean dosingInformation) {
    this.dosingInformation = dosingInformation;
  }


  public boolean isOtherPrescribingGuidance() {
    return this.otherPrescribingGuidance;
  }

  public void setOtherPrescribingGuidance(boolean otherPrescribingGuidance) {
    this.otherPrescribingGuidance = otherPrescribingGuidance;
  }


  public Markdown getText() {
    return text;
  }

  public void setText(Markdown text) {
    this.text = text;
  }


  public List<String> getImplications() {
    return implications;
  }

  public void setImplications(List<String> implications) {
    this.implications = implications;
  }


  public Map<String,Object> getLookupKey() {
    return lookupKey;
  }

  public void setLookupKey(Map<String,Object> lookupKey) {
    this.lookupKey = lookupKey;
  }


  @Override
  public int compareTo(@Nonnull RecommendationAnnotation o) {

    if (id == null) {
      return -1;
    }
    else if (o.id == null) {
      return 1;
    }
    else {
      return id.compareTo(o.getId());
    }
  }

  public boolean appliesToDrug(String drugName) {
    return relatedChemicals != null && relatedChemicals.stream()
        .anyMatch(c -> c.getName().equalsIgnoreCase(drugName));
  }

  @Override
  public String toString() {
    return getName();
  }

  /**
   * Tests whether this recommendation annotation matches the given {@link Genotype}.
   * @param genotype a possible Genotype of the subject
   * @return true if the genotype matches this recommendation
   */
  public boolean matchesGenotype(Genotype genotype) {
    if (getLookupKey() == null) {
      return false;
    }
    return genotype.getLookupKeys().stream().anyMatch(k -> RecommendationUtils.mapContains(k, getLookupKey()));
  }
}

package org.pharmgkb.pharmcat.reporter.model.result;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.stream.Collectors;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;
import org.pharmgkb.pharmcat.reporter.model.cpic.Recommendation;
import org.pharmgkb.pharmcat.reporter.model.pgkb.Group;


public class AnnotationGroup {
  @Expose
  @SerializedName("implications")
  private final Map<String,String> implications;
  @Expose
  @SerializedName("drugRecommendation")
  private final String drugRecommendation;
  @Expose
  @SerializedName("classification")
  private final String classification;
  @Expose
  @SerializedName("phenotypes")
  private final Map<String,String> phenotypes;
  @Expose
  @SerializedName("activityScore")
  private final Map<String,String> activityScore;
  @Expose
  @SerializedName("population")
  private final String population;
  @Expose
  @SerializedName("genotypes")
  private final List<Genotype> genotypes = new ArrayList<>();

  public AnnotationGroup(Recommendation recommendation) {
    implications = new TreeMap<>();
    if (recommendation.getImplications() != null) {
      implications.putAll(recommendation.getImplications());
    }
    drugRecommendation = recommendation.getDrugRecommendation();
    classification = recommendation.getClassification();
    phenotypes = new TreeMap<>();
    if (recommendation.getPhenotypes() != null) {
      phenotypes.putAll(recommendation.getPhenotypes());
    }
    activityScore = new TreeMap<>();
    if (recommendation.getActivityScore() != null) {
      activityScore.putAll(recommendation.getActivityScore());
    }
    population = recommendation.getPopulation();
  }

  public AnnotationGroup(Group group, String gene) {
    implications = new TreeMap<>();
    if (group.getImplications() != null) {
      implications.put(gene, group.getImplications().getHtmlStripped());
    }
    if (group.getRecommendation() != null) {
      drugRecommendation = group.getRecommendation().getHtmlStripped();
    } else {
      drugRecommendation = null;
    }
    if (group.getStrength() != null) {
      classification = group.getStrength().getTerm();
    } else {
      classification = null;
    }
    activityScore = new TreeMap<>();
    if (group.getActivityScore() != null) {
      activityScore.put(gene, group.getActivityScore().getHtmlStripped());
    }
    population = "";
    phenotypes = new TreeMap<>(); //TODO: fill this in?
  }

  public void addGenotype(Genotype genotype) {
    genotypes.add(genotype);
  }

  public String getDrugRecommendation() {
    return drugRecommendation;
  }

  public String getClassification() {
    return classification;
  }

  public String getPopulation() {
    return population;
  }

  public String getImplications() {
    return implications.keySet().stream()
        .map(k -> k + ": " + implications.get(k))
        .collect(Collectors.joining("\n"));
  }

  public String getPhenotypes() {
    return phenotypes.keySet().stream()
        .map(k -> k + ": " + phenotypes.get(k))
        .collect(Collectors.joining("\n"));
  }

  public String getActivityScores() {
    return activityScore.keySet().stream()
        .map(k -> k + ": " + activityScore.get(k))
        .collect(Collectors.joining("\n"));
  }

  public List<Genotype> getGenotypes() {
    return genotypes;
  }
}

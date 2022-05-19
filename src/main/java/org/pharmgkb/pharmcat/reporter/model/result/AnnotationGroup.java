package org.pharmgkb.pharmcat.reporter.model.result;

import java.sql.Array;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
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
    implications = new HashMap<>();
    if (recommendation.getImplications() != null) {
      implications.putAll(recommendation.getImplications());
    }
    drugRecommendation = recommendation.getDrugRecommendation();
    classification = recommendation.getClassification();
    phenotypes = new HashMap<>();
    if (recommendation.getPhenotypes() != null) {
      phenotypes.putAll(recommendation.getPhenotypes());
    }
    activityScore = new HashMap<>();
    if (recommendation.getActivityScore() != null) {
      activityScore.putAll(recommendation.getActivityScore());
    }
    population = recommendation.getPopulation();
  }

  public AnnotationGroup(Group group, String gene) {
    implications = new HashMap<>();
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
    activityScore = new HashMap<>();
    if (group.getActivityScore() != null) {
      activityScore.put(gene, group.getActivityScore().getHtmlStripped());
    }
    population = "";
    phenotypes = new HashMap<>(); //TODO: fill this in?
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

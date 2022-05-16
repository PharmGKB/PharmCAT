package org.pharmgkb.pharmcat.reporter.model.result;

import java.sql.Array;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.pharmgkb.pharmcat.reporter.model.cpic.Recommendation;
import org.pharmgkb.pharmcat.reporter.model.pgkb.Group;


public class AnnotationGroup {
  private final Map<String,String> implications;
  private final String drugRecommendation;
  private final String classification;
  private final Map<String,String> phenotypes;
  private final Map<String,String> activityScore;
  private final String population;
  private final List<Map<String,String>> genotypes;

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
    genotypes = new ArrayList<>();
    if (recommendation.getGenotypes() != null) {
      genotypes.addAll(recommendation.getGenotypes());
    }
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

    genotypes = new ArrayList<>();
    Map<String,String> genotype = new HashMap<>();
    genotype.put(gene, String.join("; ", group.getGenePhenotypes()));
    genotypes.add(genotype);
  }
}

package org.pharmgkb.pharmcat.reporter.model.result;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.stream.Collectors;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;
import org.pharmgkb.pharmcat.reporter.TextConstants;
import org.pharmgkb.pharmcat.reporter.model.MessageAnnotation;
import org.pharmgkb.pharmcat.reporter.model.cpic.Recommendation;
import org.pharmgkb.pharmcat.reporter.model.pgkb.Group;


public class AnnotationReport {
  @Expose
  @SerializedName("implications")
  private final SortedMap<String,String> implications = new TreeMap<>();
  @Expose
  @SerializedName("drugRecommendation")
  private String drugRecommendation;
  @Expose
  @SerializedName("classification")
  private String classification;
  @Expose
  @SerializedName("phenotypes")
  private final SortedMap<String,String> phenotypes = new TreeMap<>();
  @Expose
  @SerializedName("activityScore")
  private final SortedMap<String,String> activityScore = new TreeMap<>();
  @Expose
  @SerializedName("population")
  private String population = TextConstants.NA;
  @Expose
  @SerializedName("genotypes")
  private final List<Genotype> genotypes = new ArrayList<>();
  @Expose
  @SerializedName("comments")
  private String comments = TextConstants.NA;
  @Expose
  @SerializedName("messages")
  private final List<MessageAnnotation> m_messages = new ArrayList<>();
  @Expose
  @SerializedName("highlightedVariants")
  private final List<String> m_highlightedVariants = new ArrayList<>();




  /**
   * Create a new {@link AnnotationReport} from a CPIC {@link Recommendation}.
   *
   * @param recommendation a CPIC recommendation
   */
  public AnnotationReport(Recommendation recommendation) {
    if (recommendation.getImplications() != null) {
      implications.putAll(recommendation.getImplications());
    }
    drugRecommendation = recommendation.getDrugRecommendation();
    classification = recommendation.getClassification();
    comments = recommendation.getComments();
    if (recommendation.getPhenotypes() != null) {
      phenotypes.putAll(recommendation.getPhenotypes());
    }
    if (recommendation.getActivityScore() != null) {
      activityScore.putAll(recommendation.getActivityScore());
    }
    population = recommendation.getPopulation();
  }

  /**
   * Create a new {@link AnnotationReport} from a DPWG/PharmGKB {@link Group}.
   *
   * @param group a group of DPWG annotations
   * @param gene the single gene this annotation applies to, fix for mapping certain annotations
   */
  public AnnotationReport(Group group, String gene) {
    if (group.getImplications() != null) {
      implications.put(gene, group.getImplications().getHtmlStripped());
    }
    if (group.getRecommendation() != null) {
      drugRecommendation = group.getRecommendation().getHtmlStripped();
    }
    if (group.getStrength() != null) {
      classification = group.getStrength().getTerm();
    }
    if (group.getActivityScore() != null) {
      activityScore.put(gene, group.getActivityScore().getHtmlStripped());
    }
    if (group.getMetabolizerStatus() != null) {
      phenotypes.put(gene, group.getMetabolizerStatus().getHtmlStripped());
    }
    population = TextConstants.NA;
    comments = TextConstants.NA;
  }

  private AnnotationReport() {
  }

  /**
   * Creates a special {@link AnnotationReport} for warfarin in CPIC.
   */
  public static AnnotationReport forWarfarin(List<Genotype> genotypes) {
    AnnotationReport annotationReport = new AnnotationReport();
    genotypes.forEach(annotationReport::addGenotype);
    return annotationReport;
  }


  public List<Genotype> getGenotypes() {
    return genotypes;
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

  public Map<String, String> getPhenotypes() {
    return phenotypes;
  }

  public Map<String, String> getActivityScores() {
    return activityScore;
  }

  public String getComments() {
    return comments;
  }


  public List<MessageAnnotation> getMessages(){
    return m_messages;
  }

  public void addMessage(MessageAnnotation message) {
    m_messages.add(message);
  }


  public List<String> getHighlightedVariants() {
    return m_highlightedVariants;
  }

  public void addHighlightedVariant(String var) {
    m_highlightedVariants.add(var);
  }
}

package org.pharmgkb.pharmcat.reporter.model.result;

import java.util.List;
import java.util.Map;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;
import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;
import org.checkerframework.checker.nullness.qual.NonNull;
import org.pharmgkb.common.util.ComparisonChain;
import org.pharmgkb.pharmcat.reporter.TextConstants;
import org.pharmgkb.pharmcat.reporter.model.MessageAnnotation;
import org.pharmgkb.pharmcat.reporter.model.cpic.Recommendation;
import org.pharmgkb.pharmcat.reporter.model.pgkb.Group;


public class AnnotationReport implements Comparable<AnnotationReport> {
  @Expose
  @SerializedName("implications")
  private final SortedMap<String, String> m_implications = new TreeMap<>();
  @Expose
  @SerializedName("drugRecommendation")
  private String m_drugRecommendation;
  @Expose
  @SerializedName("classification")
  private String m_classification;
  @Expose
  @SerializedName("phenotypes")
  private final SortedMap<String, String> m_phenotypes = new TreeMap<>();
  @Expose
  @SerializedName("activityScore")
  private final SortedMap<String, String> m_activityScore = new TreeMap<>();
  @Expose
  @SerializedName("population")
  private String m_population = TextConstants.NA;
  @Expose
  @SerializedName("genotypes")
  private final SortedSet<Genotype> m_genotypes = new TreeSet<>();
  @Expose
  @SerializedName("comments")
  private String m_comments = TextConstants.NA;
  @Expose
  @SerializedName("messages")
  private final SortedSet<MessageAnnotation> m_messages = new TreeSet<>();
  @Expose
  @SerializedName("highlightedVariants")
  private final SortedSet<String> m_highlightedVariants = new TreeSet<>();


  /**
   * ID for this report.
   * Mainly used for testing.
   */
  private transient final String m_localId;

  /**
   * Gets a local ID for this annotation.
   * This should only be used to disambiguate annotations.
   */
  public String getLocalId() {
    return m_localId;
  }


  /**
   * Create a new {@link AnnotationReport} from a CPIC {@link Recommendation}.
   *
   * @param recommendation a CPIC recommendation
   */
  public AnnotationReport(Recommendation recommendation, String localId) {
    if (recommendation.getImplications() != null) {
      m_implications.putAll(recommendation.getImplications());
    }
    m_drugRecommendation = recommendation.getDrugRecommendation();
    m_classification = recommendation.getClassification();
    m_comments = recommendation.getComments();
    if (recommendation.getPhenotypes() != null) {
      m_phenotypes.putAll(recommendation.getPhenotypes());
    }

    if (recommendation.getActivityScore() != null) {
      m_activityScore.putAll(recommendation.getActivityScore());
    }
    m_population = recommendation.getPopulation();

    m_localId = localId;
  }

  /**
   * Create a new {@link AnnotationReport} from a DPWG/PharmGKB {@link Group}.
   *
   * @param group a group of DPWG annotations
   * @param gene the single gene this annotation applies to, fix for mapping certain annotations
   */
  public AnnotationReport(Group group, String gene, String localId) {
    if (group.getImplications() != null) {
      m_implications.put(gene, group.getImplications().getHtmlStripped());
    }
    if (group.getRecommendation() != null) {
      m_drugRecommendation = group.getRecommendation().getHtmlStripped();
    }
    if (group.getStrength() != null) {
      m_classification = group.getStrength().getTerm();
    }
    if (group.getActivityScore() != null) {
      m_activityScore.put(gene, group.getActivityScore().getHtmlStripped());
    }
    if (group.getMetabolizerStatus() != null) {
      m_phenotypes.put(gene, group.getMetabolizerStatus().getHtmlStripped());
    }
    m_population = TextConstants.NA;
    m_comments = TextConstants.NA;

    m_localId = localId;
  }

  private AnnotationReport(String localId) {
    m_localId = localId;
  }

  /**
   * Creates a special {@link AnnotationReport} for warfarin in CPIC.
   */
  public static AnnotationReport forWarfarin(List<Genotype> genotypes) {
    AnnotationReport annotationReport = new AnnotationReport("warfarin-cpic-1-1");
    genotypes.forEach(annotationReport::addGenotype);
    return annotationReport;
  }


  public SortedSet<Genotype> getGenotypes() {
    return m_genotypes;
  }

  public void addGenotype(Genotype genotype) {
    m_genotypes.add(genotype);
  }

  public String getDrugRecommendation() {
    return m_drugRecommendation;
  }

  public String getClassification() {
    return m_classification;
  }

  public String getPopulation() {
    return m_population;
  }

  public Map<String, String> getImplications() {
    return m_implications;
  }

  public Map<String, String> getPhenotypes() {
    return m_phenotypes;
  }

  public Map<String, String> getActivityScores() {
    return m_activityScore;
  }

  public String getComments() {
    return m_comments;
  }


  public SortedSet<MessageAnnotation> getMessages(){
    return m_messages;
  }

  public void addMessage(MessageAnnotation message) {
    m_messages.add(message);
  }


  public SortedSet<String> getHighlightedVariants() {
    return m_highlightedVariants;
  }

  public void addHighlightedVariant(String var) {
    m_highlightedVariants.add(var);
  }


  /**
   * Checks diplotypes for overriding phenotype from outside call.
   * If found, replaces phenotype from recommendation.
   */
  public void checkDiplotypes() {
    Multimap<String, String> mismatches = HashMultimap.create();
    for (Genotype genotype : m_genotypes) {
      for (Diplotype diplotype : genotype.getDiplotypes()) {
        String expected = diplotype.getOutsidePhenotypeMismatch();
        if (expected != null) {
          mismatches.put(diplotype.getGene(), diplotype.getPhenotypes().get(0));
        }
      }
    }
    for (String gene : mismatches.keySet()) {
      m_phenotypes.put(gene, String.join("/", mismatches.values()));
    }
  }


  @Override
  public int compareTo(@NonNull AnnotationReport o) {
    if (o == this) {
      return 0;
    }
    return new ComparisonChain()
        .compare(m_localId, o.getLocalId())
        .compare(m_genotypes, o.getGenotypes())
        .compare(m_drugRecommendation, o.getDrugRecommendation())
        .compare(m_classification, o.getClassification())
        .compare(m_population, o.getPopulation())
        .compare(m_comments, o.getComments())
        .compare(m_highlightedVariants, o.getHighlightedVariants())
        .compare(m_messages, o.getMessages())
        .compare(m_implications, o.getImplications())
        .compare(m_phenotypes, o.getPhenotypes())
        .compare(m_activityScore, o.getActivityScores())
        .result();
  }
}

package org.pharmgkb.pharmcat.reporter.model.result;

import java.util.ArrayList;
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
import org.pharmgkb.pharmcat.UnexpectedStateException;
import org.pharmgkb.pharmcat.reporter.TextConstants;
import org.pharmgkb.pharmcat.reporter.model.MessageAnnotation;
import org.pharmgkb.pharmcat.reporter.model.pgkb.RecommendationAnnotation;


public class AnnotationReport implements Comparable<AnnotationReport> {
  @Expose
  @SerializedName("implications")
  private final List<String> m_implications = new ArrayList<>();
  @Expose
  @SerializedName("drugRecommendation")
  private String m_drugRecommendation;
  @Expose
  @SerializedName("classification")
  private String m_classification;
  @Expose
  @SerializedName("activityScore")
  private final SortedMap<String, String> m_activityScores = new TreeMap<>();
  @Expose
  @SerializedName("population")
  private String m_population = TextConstants.NA;
  @Expose
  @SerializedName("genotypes")
  private final SortedSet<Genotype> m_genotypes = new TreeSet<>();
  @Expose
  @SerializedName("messages")
  private final SortedSet<MessageAnnotation> m_messages = new TreeSet<>();
  @Expose
  @SerializedName("highlightedVariants")
  private final SortedSet<String> m_highlightedVariants = new TreeSet<>();
  @Expose
  @SerializedName("dosingInformation")
  private boolean m_dosingInformation = false;
  @Expose
  @SerializedName("alternateDrugAvailable")
  private boolean m_alternateDrugAvailable = false;
  @Expose
  @SerializedName("otherPrescribingGuidance")
  private boolean m_otherPrescribingGuidance = false;
  @Expose
  @SerializedName("phenotypes")
  private Map<String,String> m_phenotypes = new TreeMap<>();
  @Expose
  @SerializedName("lookupKey")
  private Map<String,Object> m_lookupKey = new TreeMap<>();


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
   * Create a new {@link AnnotationReport} from a {@link RecommendationAnnotation}.
   *
   * @param recommendation a CPIC recommendation
   */
  public AnnotationReport(RecommendationAnnotation recommendation, String localId) {
    if (recommendation.getImplications() != null) {
      m_implications.addAll(recommendation.getImplications());
    }
    m_drugRecommendation = recommendation.getText().getHtmlStripped();
    m_classification = recommendation.getClassification() != null ? recommendation.getClassification().getTerm() : null;

    m_population = recommendation.getPopulation();

    m_localId = localId;

    m_dosingInformation = recommendation.isDosingInformation();
    m_alternateDrugAvailable = recommendation.isAlternateDrugAvailable();
    m_otherPrescribingGuidance = recommendation.isOtherPrescribingGuidance();
    m_lookupKey.putAll(recommendation.getLookupKey());
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

    for (Diplotype dip : genotype.getDiplotypes()) {
      String geneSymbol = dip.getGene();
      for (String phenotype : dip.getPhenotypes()) {
        String oldPhenotype = m_phenotypes.put(geneSymbol, phenotype);
        if (!dip.isAllelePresenceType() && oldPhenotype != null && !oldPhenotype.equals(phenotype)) {
          throw new UnexpectedStateException("Multiple phenotypes for gene " + geneSymbol);
        }
      }
      if (genotype.usesActivityScore()) {
        // if this genotype uses activity score at all, then we need a value for activity score for each diplotype in
        // this genotype
        if (dip.isActivityScoreType()) {
          // if the diplotype uses AS, try to use the AS in the diplotype and fallback to the AS in the recommendation
          // if it's not available
          String activityScore = TextConstants.isUnspecified(dip.getActivityScore()) ? (String)m_lookupKey.get(geneSymbol) : dip.getActivityScore();
          String oldActivity = m_activityScores.put(geneSymbol, activityScore);
          if (oldActivity != null && !oldActivity.equals(activityScore)) {
            throw new UnexpectedStateException("Multiple activity scores for gene " + geneSymbol);
          }
        } else {
          /* if this diplotype does not use AS then just mark it as n/a */
          m_activityScores.put(geneSymbol, TextConstants.NA);
        }
      }
    }
  }

  public Map<String,String> getPhenotypes() {
    return m_phenotypes;
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

  public List<String> getImplications() {
    return m_implications;
  }

  public Map<String, String> getActivityScores() {
    return m_activityScores;
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

  public boolean isDosingInformation() {
    return m_dosingInformation;
  }

  public boolean isOtherPrescribingGuidance() {
    return m_otherPrescribingGuidance;
  }

  public boolean isAlternateDrugAvailable() {
    return m_alternateDrugAvailable;
  }

  public boolean hasTags() {
    return m_dosingInformation || m_alternateDrugAvailable || m_otherPrescribingGuidance;
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
  }


  @Override
  public int compareTo(@NonNull AnnotationReport o) {
    if (o == this) {
      return 0;
    }
    return new ComparisonChain()
        .compare(m_genotypes, o.getGenotypes())
        .compare(m_population, o.getPopulation())
        .compare(m_highlightedVariants, o.getHighlightedVariants())
        .compare(m_activityScores, o.getActivityScores())
        .compare(m_classification, o.getClassification())
        .compare(m_drugRecommendation, o.getDrugRecommendation())
        .compare(m_implications, o.getImplications())
        .compare(m_messages, o.getMessages())
        .compare(m_localId, o.getLocalId())
        .result();
  }
}

package org.pharmgkb.pharmcat.reporter.model.result;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Date;
import java.util.List;
import java.util.Objects;
import java.util.Set;
import java.util.TreeSet;
import com.google.common.base.Preconditions;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableSet;
import com.google.common.collect.Multimap;
import com.google.common.collect.TreeMultimap;
import org.apache.commons.lang3.StringUtils;
import org.checkerframework.checker.nullness.qual.Nullable;
import org.pharmgkb.pharmcat.reporter.model.MessageAnnotation;
import org.pharmgkb.pharmcat.reporter.model.cpic.Drug;
import org.pharmgkb.pharmcat.reporter.model.cpic.Publication;
import org.pharmgkb.pharmcat.reporter.model.cpic.Recommendation;


/**
 * This class is a wrapper around the {@link Drug} class that also handles the matching of genotype
 * functions to recommendations.
 *
 * @author Ryan Whaley
 */
public class DrugReport implements Comparable<DrugReport> {
  private static final List<String> sf_notApplicableMatches = ImmutableList.of("PA166104949");

  private final Drug m_drug;
  private boolean m_reportable = false;
  private boolean m_isIncidentalResult = false;
  private final List<Recommendation> m_matchingRecommendations = new ArrayList<>();
  private final Multimap<String,String> m_matchedDiplotypes = TreeMultimap.create();
  private final Set<String> m_uncalledGenes = new TreeSet<>();
  private final List<MessageAnnotation> m_messages = new ArrayList<>();
  private final List<String> m_reportVariants = new ArrayList<>();

  public DrugReport(Drug drug) {
    m_drug = drug;
  }

  /**
   * Gets the name of the guideline, a pass-through to the stored guideline.
   */
  public String getName() {
    return m_drug.getGuidelineName();
  }

  /**
   * Gets the ID of the guideline
   */
  public String getId() {
    return m_drug.getDrugId();
  }

  /**
   * Gets just the symbols of the related genes of the guideline. Calculated from data in the original guideline.
   */
  public List<String> getRelatedGeneSymbols() {
    return m_drug.getGenes();
  }

  public Set<String> getRelatedDrugs() {
    return ImmutableSet.of(m_drug.getDrugName());
  }

  public List<Recommendation> getRecommendations() {
    if (m_drug == null || m_drug.getRecommendations() == null) {
      return new ArrayList<>();
    } else {
      return m_drug.getRecommendations();
    }
  }

  public boolean isMatched() {
    return sf_notApplicableMatches.contains(getId()) || m_matchingRecommendations.size()>0;
  }

  public boolean hasMultipleMatches() {
    return m_matchingRecommendations.size()>1;
  }


  public List<Recommendation> getMatchingRecommendations() {
    return m_matchingRecommendations;
  }

  public void addMatchingRecommendation(Recommendation recommendation) {
    m_matchingRecommendations.add(recommendation);
  }

  /**
   * Gets the URL for the whole annotation
   */
  public String getUrl() {
    return m_drug.getUrl();
  }

  /**
   * Gets the matched diplotypes for this annotation (should be a subset of all called genotypes).
   *
   * Will be in functional form, e.g. "GENEX:No Function/Increased Function"
   */
  public Multimap<String,String> getMatchedDiplotypes() {
    return m_matchedDiplotypes;
  }

  /**
   * True if each of the genes in this guideline has at least one called diplotype
   */
  public boolean isReportable() {
    return sf_notApplicableMatches.contains(getId()) || m_reportable;
  }

  public void setReportable(boolean reportable) {
    m_reportable = reportable;
  }

  @Override
  public int compareTo(DrugReport o) {
    int rez = Boolean.compare(isReportable(), o.isReportable());
    if (rez != 0) {
      return rez * -1;
    }
    rez = Objects.compare(getName(), o.getName(), String.CASE_INSENSITIVE_ORDER);
    return rez;
  }

  public Set<String> getUncalledGenes() {
    return m_uncalledGenes;
  }

  public void addUncalledGene(String geneSymbol) {
    m_uncalledGenes.add(geneSymbol);
  }

  /**
   * TODO(ryan): needs to be re-implemented
   * @return true if this guideline matches to a recommendation with an Rx Change
   */
  public boolean isRxChange() {
    return false;
  }

  /**
   * TODO(ryan): needs to be re-implemented
   * @return true if this guideline matches to a recommendation with a "Possibly" Rx Change
   */
  public boolean isRxPossible() {
    return false;
  }

  /**
   * Finds the matching {@link Recommendation} objects for the given <code>reportGenotype</code>, adds it to the group, and then
   * marks it as a match.
   * @param reportGenotype a multi-gene genotype function String in the form of "GENEA:No Function/No Function;GENEB:Normal Function/Normal Function"
   */
  public void addReportGenotype(String reportGenotype) {
    Preconditions.checkArgument(StringUtils.isNotBlank(reportGenotype));

    getRecommendations().stream()
        .filter(r -> r.matchLookupKey(reportGenotype))
        .forEach(r -> {
          addMatchingRecommendation(r);
          r.addMatchedDiplotype(reportGenotype);
        });
  }

  public List<MessageAnnotation> getMessages() {
    return m_messages;
  }

  public void addMessages(@Nullable Collection<MessageAnnotation> messages) {
    if (messages == null) {
      return;
    }

    // separate the general messages from specific genotype call messages
    messages.forEach(ma -> {
      if (ma.getExceptionType().equals(MessageAnnotation.TYPE_GENOTYPE)) {
        m_reportVariants.add(ma.getMatches().getVariant());
      }
      else {
        m_messages.add(ma);
      }
    });
  }

  /**
   * Gets whether any related gene has any incidental allele called for it
   */
  public boolean isIncidentalResult() {
    return m_isIncidentalResult;
  }

  public void setIncidentalResult(boolean incidentalResult) {
    m_isIncidentalResult = incidentalResult;
  }

  public List<String> getReportVariants() {
    return m_reportVariants;
  }

  public String toString() {
    return String.join(", ", getRelatedDrugs());
  }

  /**
   * Gets the literature objects that are used for citation of this Guideline
   */
  public List<Publication> getCitations() {
    return m_drug.getPublications();
  }

  /**
   * Date of last modification for the guideline information.
   *
   * This is currently not supported by CPIC-sourced guideline data, will eventually re-implement.
   *
   * @return a nullable Date of last modification
   */
  public Date getLastModified() {
    return null;
  }
}

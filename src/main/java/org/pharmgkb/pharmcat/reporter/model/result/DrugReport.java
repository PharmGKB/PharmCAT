package org.pharmgkb.pharmcat.reporter.model.result;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Date;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import com.google.common.base.MoreObjects;
import com.google.common.base.Preconditions;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableSet;
import org.checkerframework.checker.nullness.qual.Nullable;
import org.pharmgkb.pharmcat.reporter.model.DataSource;
import org.pharmgkb.pharmcat.reporter.model.MessageAnnotation;
import org.pharmgkb.pharmcat.reporter.model.cpic.Drug;
import org.pharmgkb.pharmcat.reporter.model.cpic.Publication;
import org.pharmgkb.pharmcat.reporter.model.cpic.Recommendation;
import org.pharmgkb.pharmcat.reporter.model.pgkb.GuidelinePackage;


/**
 * This class is a wrapper around the {@link Drug} class that also handles the matching of genotype
 * functions to recommendations.
 *
 * @author Ryan Whaley
 */
public class DrugReport implements Comparable<DrugReport> {
  // NOTE: This is so the "No Recommendations" section doesn't show in the warfarin guideline
  private static final List<String> sf_notApplicableMatches = ImmutableList.of("RxNorm:11289"); // ID for warfarin

  private final String f_name;
  private String m_cpicId;
  private String m_pgkbId;
  private String m_cpicVersion;
  private boolean m_reportable = false;
  private final Set<String> f_relatedGenes = new TreeSet<>();
  private final List<Recommendation> f_allRecommendations = new ArrayList<>();
  private final Set<Recommendation> m_matchingRecommendations = new HashSet<>();
  private final Set<String> m_uncalledGenes = new TreeSet<>();
  private final List<MessageAnnotation> m_messages = new ArrayList<>();
  private final List<String> m_reportVariants = new ArrayList<>();
  private final List<String> f_urls = new ArrayList<>();
  private final List<Publication> f_publications = new ArrayList<>();
  private DataSource m_source;
  private List<GuidelinePackage> f_dpwgGuidelinePackages = new ArrayList<>();
  private final List<GuidelineReport> f_guidelines = new ArrayList<>();

  public DrugReport(String name) {
    f_name = name;
  }

  public GuidelineReport addDrugData(Drug drug) {
    m_cpicVersion = drug.getCpicVersion();
    m_cpicId = drug.getDrugId();
    Preconditions.checkArgument(f_name.equalsIgnoreCase(drug.getDrugName()));
    if (drug.getRecommendations() != null) {
      f_allRecommendations.addAll(drug.getRecommendations());
    }
    f_relatedGenes.addAll(drug.getGenes());
    f_urls.add(drug.getUrl());
    if (drug.getPublications() != null) {
      f_publications.addAll(drug.getPublications());
    }
    m_source = DataSource.CPIC;

    GuidelineReport guidelineReport = new GuidelineReport(drug);
    addGuideline(guidelineReport);
    return guidelineReport;
  }

  public GuidelineReport addDrugData(GuidelinePackage guidelinePackage) {
    addDpwgGuidelinePackage(guidelinePackage);
    m_pgkbId = guidelinePackage.getGuideline().getId();
    guidelinePackage.getGuideline().getGuidelineGenes().stream()
        .map(gg -> gg.getGene().getSymbol())
        .forEach(f_relatedGenes::add);
    f_urls.add(guidelinePackage.getGuideline().getUrl());
    m_source = DataSource.DPWG;

    GuidelineReport guidelineReport = new GuidelineReport(guidelinePackage);
    guidelinePackage.getMatchedGroups()
        .forEach((group) -> guidelineReport.addAnnotationGroup(new AnnotationGroup(group, guidelinePackage.getGenes().iterator().next())));
    addGuideline(guidelineReport);
    return guidelineReport;
  }

  /**
   * Gets the name of the guideline, a pass-through to the stored guideline.
   */
  public String getName() {
    return f_name;
  }

  /**
   * Gets the version of CPIC the resource data is based on
   * @return a String describing the CPIC version
   */
  public String getCpicVersion() {
    return m_cpicVersion;
  }

  /**
   * Gets the CPIC or the PGKB ID, whichever is specified, in that order
   * @return an ID for this drug
   */
  public String getId() {
    return MoreObjects.firstNonNull(getCpicId(), getPgkbId());
  }

  /**
   * Gets the CPIC ID of the drug
   */
  public String getCpicId() {
    return m_cpicId;
  }

  /**
   * Gets the PharmGKB ID of the drug
   */
  public String getPgkbId() {
    return m_pgkbId;
  }

  /**
   * Gets just the symbols of the related genes of the guideline. Calculated from data in the original guideline.
   */
  public Collection<String> getRelatedGeneSymbols() {
    return f_relatedGenes;
  }

  public Set<String> getRelatedDrugs() {
    return ImmutableSet.of(f_name);
  }

  public List<Recommendation> getRecommendations() {
    return f_allRecommendations;
  }

  public boolean isMatched() {
    return sf_notApplicableMatches.contains(getCpicId())
        || m_matchingRecommendations.size()>0
        || f_dpwgGuidelinePackages.stream().anyMatch(GuidelinePackage::hasMatch);
  }

  /**
   * Determines if a DrugReport is "ignored". If ignored, it will not be displayed in the reporter output. The drug
   * report is ignored if any gene used for the recommendations is ignored as determined by
   * {@link GeneReport#isIgnored(String)}.
   * @return true if the drug report will not be included
   */
  public boolean isIgnored() {
    return f_relatedGenes.stream().anyMatch(GeneReport::isIgnored);
  }

  public boolean hasMultipleMatches() {
    return m_matchingRecommendations.size()>1;
  }


  public Set<Recommendation> getMatchingRecommendations() {
    return m_matchingRecommendations;
  }

  private void addMatchingRecommendation(Recommendation recommendation) {
    m_matchingRecommendations.add(recommendation);
  }

  /**
   * Gets the URL for the whole annotation
   */
  public List<String> getUrls() {
    return f_urls;
  }

  /**
   * True if each of the genes in this guideline has at least one called diplotype
   */
  public boolean isReportable() {
    return sf_notApplicableMatches.contains(getCpicId()) || m_reportable;
  }

  public void setReportable(boolean reportable) {
    m_reportable = reportable;
  }

  @Override
  public int compareTo(DrugReport o) {
    int rez = toString().compareToIgnoreCase(o.toString());
    if (rez != 0) {
      return rez;
    }
    rez = Boolean.compare(isReportable(), o.isReportable());
    return rez;
  }

  public Set<String> getUncalledGenes() {
    return m_uncalledGenes;
  }

  public void addUncalledGene(String geneSymbol) {
    m_uncalledGenes.add(geneSymbol);
  }

  /**
   * Finds the matching {@link Recommendation} objects for the given <code>phenotypeKey</code>, adds it to the group,
   * and then marks it as a match.
   * @param phenotypeKey a Map of gene symbol to phenotype String
   * @param diplotypes a Set of diplotypes strings
   */
  public void addReportGenotype(Map<String,String> phenotypeKey, SortedSet<String> diplotypes) {
    Preconditions.checkNotNull(phenotypeKey);
    Preconditions.checkArgument(!phenotypeKey.isEmpty());

    getRecommendations().stream()
        .filter((r) -> r.getLookupKey().equals(phenotypeKey))
        .forEach((r) -> {
          addMatchingRecommendation(r);
          r.addMatchedDiplotype(String.join(", ", diplotypes));
          f_guidelines.stream()
              .filter((g) -> g.getSource() == DataSource.CPIC)
              .forEach((g) -> {
                g.addAnnotationGroup(new AnnotationGroup(r));
              });
        });
  }

  /**
   * Experimental new genotype matcher
   * @param genotype
   */
  public void matchAnnotationsToGenotype(Genotype genotype) {
    getRecommendations().stream()
        .filter((r) -> r.matchesGenotype(genotype))
        .forEach((r) -> {
          addMatchingRecommendation(r);
          r.addMatchedGenotype(genotype);
          f_guidelines.stream()
              .filter((g) -> g.getSource() == DataSource.CPIC)
              .forEach((g) -> {
                g.addAnnotationGroup(new AnnotationGroup(r));
              });
        });
  }

  public List<MessageAnnotation> getMessages() {
    return m_messages;
  }

  public void addMessage(MessageAnnotation message) {
    m_messages.add(message);
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

  public List<String> getReportVariants() {
    return m_reportVariants;
  }

  public String toString() {
    return String.join(", ", getName());
  }

  /**
   * Gets the literature objects that are used for citation of this Guideline
   */
  public List<Publication> getCitations() {
    return f_publications;
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

  public DataSource getSource() {
    return m_source;
  }

  public void setSource(DataSource source) {
    m_source = source;
  }

  public List<GuidelinePackage> getDpwgGuidelinePackages() {
    return f_dpwgGuidelinePackages;
  }

  public void setDpwgGuidelinePackages(List<GuidelinePackage> f_dpwgGuidelinePackages) {
    this.f_dpwgGuidelinePackages = f_dpwgGuidelinePackages;
  }

  public void addDpwgGuidelinePackage(GuidelinePackage guidelinePackage) {
    f_dpwgGuidelinePackages.add(guidelinePackage);
  }

  public boolean isCpic() {
    return m_cpicId != null;
  }

  public boolean isDpwg() {
    return m_pgkbId != null;
  }

  public List<GuidelineReport> getGuidelines() {
    return f_guidelines;
  }

  public void addGuideline(GuidelineReport guidelineReport) {
    f_guidelines.add(guidelineReport);
  }
}

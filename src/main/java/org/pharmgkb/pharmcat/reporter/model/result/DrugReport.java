package org.pharmgkb.pharmcat.reporter.model.result;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;
import com.google.common.base.MoreObjects;
import com.google.common.base.Preconditions;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableSet;
import org.checkerframework.checker.nullness.qual.Nullable;
import org.pharmgkb.pharmcat.reporter.ReportContext;
import org.pharmgkb.pharmcat.reporter.model.DataSource;
import org.pharmgkb.pharmcat.reporter.model.MessageAnnotation;
import org.pharmgkb.pharmcat.reporter.model.cpic.Drug;
import org.pharmgkb.pharmcat.reporter.model.cpic.Publication;
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
  private final List<MessageAnnotation> m_messages = new ArrayList<>();
  private final List<String> m_reportVariants = new ArrayList<>();
  private final List<String> f_urls = new ArrayList<>();
  private final List<Publication> f_publications = new ArrayList<>();
  private final List<GuidelineReport> f_guidelines = new ArrayList<>();

  public DrugReport(String name) {
    f_name = name;
  }

  public void addDrugData(Drug drug, ReportContext reportContext) {
    m_cpicId = drug.getDrugId();
    Preconditions.checkArgument(f_name.equalsIgnoreCase(drug.getDrugName()));
    f_urls.add(drug.getUrl());
    if (drug.getPublications() != null) {
      f_publications.addAll(drug.getPublications());
    }

    GuidelineReport guidelineReport = new GuidelineReport(drug);
    drug.getGenes().forEach((geneSymbol) -> guidelineReport.addRelatedGeneReport(reportContext.getGeneReport(geneSymbol)));
    addGuideline(guidelineReport);
  }

  public void addDrugData(GuidelinePackage guidelinePackage, ReportContext reportContext) {
    m_pgkbId = guidelinePackage.getGuideline().getId();
    f_urls.add(guidelinePackage.getGuideline().getUrl());

    GuidelineReport guidelineReport = new GuidelineReport(guidelinePackage);
    guidelinePackage.getGenes().forEach((geneSymbol) -> guidelineReport.addRelatedGeneReport(reportContext.getGeneReport(geneSymbol)));
    guidelinePackage.getMatchedGroups()
        .forEach((group) -> guidelineReport.addAnnotationGroup(new AnnotationGroup(group, guidelinePackage.getGenes().iterator().next())));
    addGuideline(guidelineReport);
  }

  /**
   * Gets the name of the guideline, a pass-through to the stored guideline.
   */
  public String getName() {
    return f_name;
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
  private String getCpicId() {
    return m_cpicId;
  }

  /**
   * Gets the PharmGKB ID of the drug
   */
  private String getPgkbId() {
    return m_pgkbId;
  }

  /**
   * Gets just the symbols of the related genes of the guideline. Calculated from data in the original guideline.
   */
  public Collection<String> getRelatedGeneSymbols() {
    return f_guidelines.stream()
        .flatMap((guidelineReport) -> guidelineReport.getRelatedGeneReports().stream())
        .map(GeneReport::getGene)
        .distinct()
        .sorted()
        .collect(Collectors.toList());
  }

  public Set<String> getRelatedDrugs() {
    return ImmutableSet.of(f_name);
  }

  public boolean isMatched() {
    return sf_notApplicableMatches.contains(getCpicId())
        || f_guidelines.stream().anyMatch(GuidelineReport::isMatched);
  }

  /**
   * Gets the URL for the whole annotation
   */
  public List<String> getUrls() {
    return f_urls;
  }

  @Override
  public int compareTo(DrugReport o) {
    return toString().compareToIgnoreCase(o.toString());
  }

  /**
   * Experimental new genotype matcher
   * TODO: fill out docs here
   * @param genotype a possible genotype for this report
   */
  public void matchAnnotationsToGenotype(Genotype genotype, Drug cpicDrug) {
    if (cpicDrug.getRecommendations() != null) {
      cpicDrug.getRecommendations().stream()
          .filter((r) -> r.matchesGenotype(genotype))
          .forEach((r) -> f_guidelines.stream()
              .filter((g) -> g.getSource() == DataSource.CPIC)
              .forEach((g) -> {
                AnnotationGroup annGroup = new AnnotationGroup(r);
                annGroup.addGenotype(genotype);
                g.addAnnotationGroup(annGroup);
              }));
    }
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

  public long getMatchedGroupCount() {
    return getGuidelines().stream().mapToLong(g -> g.getAnnotationGroups().size()).sum();
  }
}

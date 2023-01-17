package org.pharmgkb.pharmcat.reporter.model.result;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Objects;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.stream.Collectors;
import com.google.common.base.Preconditions;
import com.google.common.collect.ImmutableList;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;
import org.checkerframework.checker.nullness.qual.NonNull;
import org.checkerframework.checker.nullness.qual.Nullable;
import org.pharmgkb.common.util.ComparisonChain;
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

  @Expose
  @SerializedName("name")
  private final String m_drugName;
  @Expose
  @SerializedName("id")
  private String m_id;
  @Expose
  @SerializedName("source")
  private DataSource m_source;
  @Expose
  @SerializedName("version")
  private String m_version;

  @Expose
  @SerializedName("messages")
  private final SortedSet<MessageAnnotation> m_messages = new TreeSet<>();
  @Expose
  @SerializedName("variants")
  private final SortedSet<String> m_reportVariants = new TreeSet<>();
  @Expose
  @SerializedName("urls")
  private final List<String> m_urls = new ArrayList<>();
  @Expose
  @SerializedName("citations")
  private final List<Publication> m_citations = new ArrayList<>();
  @Expose
  @SerializedName("guidelines")
  private final SortedSet<GuidelineReport> m_guidelines = new TreeSet<>();


  public DrugReport(Drug drug, ReportContext reportContext) {
    m_drugName = drug.getDrugName();
    m_id = drug.getDrugId();
    m_source = DataSource.CPIC;
    m_version = drug.getCpicVersion();
    m_urls.add(drug.getUrl());
    if (drug.getCitations() != null) {
      // cpic data can have array with null value in it
      drug.getCitations().stream()
          .filter(Objects::nonNull)
          .forEach(m_citations::add);
    }

    // 1 guideline report per CPIC drug
    GuidelineReport guidelineReport = new GuidelineReport(drug, reportContext);
    // link gene report back to drug report
    guidelineReport.getGeneReports().forEach((gr) -> gr.addRelatedDrug(this));
    m_guidelines.add(guidelineReport);
  }

  public DrugReport(String name, SortedSet<GuidelinePackage> guidelinePackages, ReportContext reportContext) {
    Preconditions.checkArgument(guidelinePackages != null && guidelinePackages.size() > 0);
    m_drugName = name;
    m_id = guidelinePackages.first().getGuideline().getRelatedChemicals().stream()
        .filter((c) -> c.getName().equals(name))
        .findFirst()
        .orElseThrow(() -> new IllegalStateException("DPWG guideline " +
            guidelinePackages.first().getGuideline().getId() + " is supposd to be related to " + name + " but is not"))
        .getId();
    m_source = DataSource.DPWG;
    m_version = guidelinePackages.stream().map(GuidelinePackage::getVersion).collect(Collectors.joining(", "));

    // DPWG drug can have multiple guideline reports
    for (GuidelinePackage guidelinePackage : guidelinePackages) {
      m_urls.add(guidelinePackage.getGuideline().getUrl());
      if (guidelinePackage.getCitations() != null) {
        m_citations.addAll(guidelinePackage.getCitations());
      }

      GuidelineReport guidelineReport = new GuidelineReport(guidelinePackage, reportContext);
      // link gene report back to drug report
      guidelineReport.getGeneReports().forEach((gr) -> gr.addRelatedDrug(this));
      m_guidelines.add(guidelineReport);
    }
  }


  /**
   * Gets the name of the drug this {@link DrugReport} is on.
   */
  public String getName() {
    return m_drugName;
  }

  /**
   * Gets the drug ID.
   */
  public String getId() {
    return m_id;
  }

  public DataSource getSource() {
    return m_source;
  }

  public String getVersion() {
    return m_version;
  }


  /**
   * Gets just the symbols of the related genes of the guideline. Calculated from data in the original guideline.
   */
  public Collection<String> getRelatedGeneSymbols() {
    return m_guidelines.stream()
        .flatMap((guidelineReport) -> guidelineReport.getGeneReports().stream())
        .map(GeneReport::getGene)
        .distinct()
        .sorted()
        .collect(Collectors.toList());
  }

  public boolean isMatched() {
    return (m_source == DataSource.CPIC && sf_notApplicableMatches.contains(m_id))
        || m_guidelines.stream().anyMatch(GuidelineReport::isMatched);
  }

  /**
   * Gets the URL for the whole annotation
   */
  public List<String> getUrls() {
    return m_urls;
  }

  @Override
  public int compareTo(@NonNull DrugReport o) {
    if (this == o) {
      return 0;
    }
    return new ComparisonChain()
        .compare(m_drugName, o.getName())
        .compare(m_source, o.getSource())
        .compare(m_version, o.getVersion())
        .compare(m_id, o.getId())
        .compare(m_guidelines, o.getGuidelines())
        .compare(m_messages, o.getMessages())
        .result();
  }


  public SortedSet<MessageAnnotation> getMessages() {
    return m_messages;
  }

  public boolean addMessage(MessageAnnotation message) {
    return m_messages.add(message);
  }

  public void addMessages(@Nullable Collection<MessageAnnotation> messages) {
    if (messages == null) {
      return;
    }

    // separate the general messages from specific genotype call messages
    messages.forEach(ma -> {
      if (ma.getExceptionType().equals(MessageAnnotation.TYPE_REPORT_AS_GENOTYPE)) {
        m_reportVariants.add(ma.getMatches().getVariant());
      }
      else {
        m_messages.add(ma);
      }
    });
  }

  public boolean removeMessage(MessageAnnotation message) {
    return m_messages.remove(message);
  }

  /**
   * Gets list of variants to display as part of genotype for recommendation.
   */
  public SortedSet<String> getReportVariants() {
    return m_reportVariants;
  }

  public String toString() {
    return m_drugName;
  }

  /**
   * Gets the literature objects that are used for citation of this Guideline
   */
  public List<Publication> getCitations() {
    return m_citations;
  }

  public SortedSet<GuidelineReport> getGuidelines() {
    return m_guidelines;
  }

  public int getMatchedAnnotationCount() {
    return getGuidelines().stream().mapToInt(g -> g.getAnnotations().size()).sum();
  }
}

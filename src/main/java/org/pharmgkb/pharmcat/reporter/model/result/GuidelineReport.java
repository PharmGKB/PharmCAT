package org.pharmgkb.pharmcat.reporter.model.result;

import java.util.Collection;
import java.util.List;
import java.util.Objects;
import java.util.SortedSet;
import java.util.TreeSet;
import com.google.common.collect.HashMultimap;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;
import org.checkerframework.checker.nullness.qual.NonNull;
import org.pharmgkb.common.util.ComparisonChain;
import org.pharmgkb.pharmcat.reporter.ReportContext;
import org.pharmgkb.pharmcat.reporter.model.DataSource;
import org.pharmgkb.pharmcat.reporter.model.pgkb.GuidelinePackage;
import org.pharmgkb.pharmcat.reporter.model.pgkb.RecommendationAnnotation;


public class GuidelineReport implements Comparable<GuidelineReport> {
  @Expose
  @SerializedName("id")
  private String m_id;
  @Expose
  @SerializedName("name")
  private String m_name;

  @Expose
  @SerializedName("source")
  private DataSource m_source;
  @Expose
  @SerializedName("version")
  private String m_version;
  @Expose
  @SerializedName("url")
  private String m_url;
  @Expose
  @SerializedName("annotations")
  private SortedSet<AnnotationReport> m_annotationReports = new TreeSet<>();

  private transient final SortedSet<String> m_genes = new TreeSet<>();
  private transient final SortedSet<GeneReport> m_geneReports = new TreeSet<>();
  private transient final SortedSet<Diplotype> m_sourceDiplotypes = new TreeSet<>();
  private transient List<Genotype> m_recommendationGenotypes;
  private transient final SortedSet<String> m_homozygousComponentHaplotypes = new TreeSet<>();


  /**
   * Private constructor for GSON.
   */
  @SuppressWarnings("unused")
  private GuidelineReport() {
  }

  public GuidelineReport(GuidelinePackage guidelinePackage, ReportContext reportContext, String drugName) {
    m_id = guidelinePackage.getGuideline().getId();
    m_name = guidelinePackage.getGuideline().getName();
    m_source = DataSource.valueOf(guidelinePackage.getGuideline().getSource());
    m_version = guidelinePackage.getVersion();
    m_url = guidelinePackage.getGuideline().getUrl();
    initializeGenes(guidelinePackage.getGenes(), reportContext);
    matchAnnotations(guidelinePackage, drugName);
  }

  private void initializeGenes(Collection<String> genes, ReportContext reportContext) {
    // link guideline report to gene report
    for (String geneSymbol : genes) {
      GeneReport geneReport = reportContext.getGeneReport(m_source, geneSymbol);
      if (geneReport == null) {
        continue;
      }
      m_genes.add(geneSymbol);
      m_geneReports.add(geneReport);
    }
    m_geneReports.stream()
        .flatMap(gr -> gr.getSourceDiplotypes().stream())
        .forEach(m_sourceDiplotypes::add);
    m_recommendationGenotypes = Genotype.makeGenotypes(m_geneReports);
    m_geneReports.forEach(gr -> {
      m_homozygousComponentHaplotypes.addAll(gr.getMatcherHomozygousComponentHaplotypes());
    });
  }



  public String getId() {
    return m_id;
  }

  public String getName() {
    return m_name;
  }

  public DataSource getSource() {
    return m_source;
  }

  public String getVersion() {
    return m_version;
  }

  public String getUrl() {
    return m_url;
  }


  public SortedSet<AnnotationReport> getAnnotations() {
    return m_annotationReports;
  }

  public boolean isMatched() {
    return !m_annotationReports.isEmpty();
  }


  public SortedSet<String> getGenes() {
    return m_genes;
  }

  public SortedSet<GeneReport> getGeneReports() {
    return m_geneReports;
  }

  /**
   * Checks if any {@link GeneReport} is reportable.
   */
  public boolean isReportable() {
    return m_geneReports.stream()
        .anyMatch(GeneReport::isReportable);
  }

  public SortedSet<Diplotype> getSourceDiplotypes() {
    return m_sourceDiplotypes;
  }


  public SortedSet<String> getHomozygousComponentHaplotypes() {
    return m_homozygousComponentHaplotypes;
  }


  private void matchAnnotations(GuidelinePackage guidelinePackage, String drugName) {
    HashMultimap<RecommendationAnnotation, Genotype> matchedGenotypes = HashMultimap.create();
    for (Genotype genotype : m_recommendationGenotypes) {
      guidelinePackage.getRecommendations().stream()
          .filter(Objects::nonNull)
          .filter(rec -> rec.appliesToDrug(drugName))
          .filter(rec -> rec.matchesGenotype(genotype))
          .forEach(rec -> matchedGenotypes.put(rec, genotype));
    }
    if (drugName.equals("warfarin") && m_source == DataSource.CPIC) {
      AnnotationReport ann = AnnotationReport.forWarfarin(m_recommendationGenotypes);
      m_recommendationGenotypes.forEach(ann::addGenotype);
      m_annotationReports.add(ann);
    }
    for (RecommendationAnnotation recommendationAnnotation : matchedGenotypes.keys()) {
      String id = guidelinePackage.getGuideline().getSource() + "-" + recommendationAnnotation.getId();
      AnnotationReport annotationReport = new AnnotationReport(recommendationAnnotation, id);
      matchedGenotypes.get(recommendationAnnotation).forEach(annotationReport::addGenotype);
      annotationReport.checkDiplotypes();
      m_annotationReports.add(annotationReport);
    }
  }

  @Override
  public int compareTo(@NonNull GuidelineReport o) {
    if (this == o) {
      return 0;
    }
    return new ComparisonChain()
        .compare(m_name, o.getName())
        .compare(m_source, o.getSource())
        .compare(m_version, o.getVersion())
        .compare(m_id, o.getId())
        .compare(m_annotationReports, o.getAnnotations())
        .result();
  }
}

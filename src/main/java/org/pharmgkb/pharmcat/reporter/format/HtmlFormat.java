package org.pharmgkb.pharmcat.reporter.format;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.text.SimpleDateFormat;
import java.util.*;
import java.util.stream.Collectors;
import com.github.jknack.handlebars.Handlebars;
import com.github.jknack.handlebars.helper.StringHelpers;
import com.github.jknack.handlebars.io.ClassPathTemplateLoader;
import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.google.common.collect.SortedSetMultimap;
import com.google.common.collect.TreeMultimap;
import org.apache.commons.lang3.StringUtils;
import org.pharmgkb.pharmcat.reporter.ReportContext;
import org.pharmgkb.pharmcat.reporter.handlebars.ReportHelpers;
import org.pharmgkb.pharmcat.reporter.model.DataSource;
import org.pharmgkb.pharmcat.reporter.model.MessageAnnotation;
import org.pharmgkb.pharmcat.reporter.model.VariantReport;
import org.pharmgkb.pharmcat.reporter.model.cpic.Publication;
import org.pharmgkb.pharmcat.reporter.model.result.AnnotationReport;
import org.pharmgkb.pharmcat.reporter.model.result.CallSource;
import org.pharmgkb.pharmcat.reporter.model.result.Diplotype;
import org.pharmgkb.pharmcat.reporter.model.result.DrugLink;
import org.pharmgkb.pharmcat.reporter.model.result.DrugReport;
import org.pharmgkb.pharmcat.reporter.model.result.GeneReport;
import org.pharmgkb.pharmcat.reporter.model.result.Genotype;
import org.pharmgkb.pharmcat.reporter.model.result.GuidelineReport;

import static org.pharmgkb.pharmcat.reporter.DpydCaller.isDpyd;


/**
 * An HTML-formatted version of {@link ReportContext} data.
 */
public class HtmlFormat extends AbstractFormat {
  private static final String sf_templatePrefix = "/org/pharmgkb/pharmcat/reporter";
  private static final String FINAL_REPORT = "report";
  private List<DataSource> m_sources = Lists.newArrayList(DataSource.CPIC, DataSource.DPWG);
  private boolean m_compact;
  private final boolean f_testMode;

  /**
   * Constructor. Takes the path to write to and whether this output is "test mode" or not
   * @param outputPath the path to write the data to
   * @param testMode if true will leave out some time-related metadata that may give false-positives while diffing
   * output
   */
  public HtmlFormat(Path outputPath, boolean testMode) {
    super(outputPath);
    f_testMode = testMode;
  }

  public HtmlFormat sources(List<DataSource> sources) {
    if (sources != null) {
      m_sources = sources;
    }
    return this;
  }

  public HtmlFormat compact(boolean compact) {
    m_compact = compact;
    return this;
  }


  public void write(ReportContext reportContext) throws IOException {
    Map<String, Object> reportData = compile(reportContext);

    Handlebars handlebars = new Handlebars(new ClassPathTemplateLoader(sf_templatePrefix));
    StringHelpers.register(handlebars);
    handlebars.registerHelpers(ReportHelpers.class);

    try (BufferedWriter writer = Files.newBufferedWriter(getOutputPath(), StandardCharsets.UTF_8)) {
      writer.write(handlebars.compile(FINAL_REPORT).apply(reportData));
    }
  }

  /**
   * Make a Map of data that will be used in the final report. This map will be serialized and then applied to the
   * Handlebars template.
   * @return a Map of data to serialize into JSON
   */
  private Map<String,Object> compile(ReportContext reportContext) {

    Map<String,Object> result = new HashMap<>();
    if (!f_testMode) {
      result.put("generatedOn", new SimpleDateFormat("MMMMM dd, yyyy").format(reportContext.getGeneratedOn()));
      result.put("pharmcatVersion", reportContext.getPharmcatVersion());
      result.put("cpicVersion", reportContext.getCpicVersion());
    }
    result.put("compact", m_compact);

    if (StringUtils.isNotBlank(reportContext.getTitle())) {
      result.put("title", reportContext.getTitle());
    }


    // Section I: Genotype Summary
    Set<String> totalGenes = new HashSet<>();
    Set<String> calledGenes = new HashSet<>();
    boolean hasCombo = false;
    boolean hasMessages = false;
    boolean hasMissingVariants = false;
    boolean hasUnphasedNote = false;

    SortedSetMultimap<String, GeneReport> geneReportMap = TreeMultimap.create();
    for (DataSource source : reportContext.getGeneReports().keySet()) {
      if (!m_sources.contains(source)) {
        continue;
      }
      for (GeneReport geneReport : reportContext.getGeneReports().get(source).values()) {
        // skip any genes on the blacklist
        if (geneReport.isIgnored()) {
          continue;
        }
        String symbol = geneReport.getGeneDisplay();
        totalGenes.add(symbol);
        if (!m_compact) {
          geneReportMap.put(symbol, geneReport);
        }

        // skip any uncalled genes
        boolean allVariantsMissing = geneReport.getVariantReports().stream().allMatch(VariantReport::isMissing);
        if ((!geneReport.isCalled() || allVariantsMissing) && (geneReport.getRecommendationDiplotypes().isEmpty())) {
          continue;
        }
        if (geneReport.getRelatedDrugs().size() == 0) {
          continue;
        }

        if (geneReport.isReportable()) {
          calledGenes.add(symbol);
          if (m_compact) {
            geneReportMap.put(symbol, geneReport);
          }

          if (geneReport.isMissingVariants()) {
            hasMissingVariants = true;
          }
          hasCombo = hasCombo || geneReport.getMessages().stream()
              .anyMatch(m -> m.getExceptionType().equals(MessageAnnotation.TYPE_COMBO));
          hasMessages = hasMessages || hasMessages(geneReport);
          hasUnphasedNote = hasUnphasedNote || showUnphasedNote(geneReport);
        }
      }
    }

    SortedSet<String> genes = new TreeSet<>(geneReportMap.keySet());
    List<Map<String, Object>> summaries = new ArrayList<>();
    List<GeneReport> geneReports = new ArrayList<>();
    for (String symbol : genes) {
      SortedSet<GeneReport> reports = geneReportMap.get(symbol);
      if (reports.size() > 2) {
        throw new IllegalStateException("More than 2 gene reports for " + symbol);
      }
      Optional<Map<String, Object>> opt = buildGenotypeSummary(symbol, reports);
      opt.ifPresent(summaries::add);
      if (!m_compact || opt.isPresent()) {
        geneReports.add(reports.first());
      }
    }
    result.put("genes", genes);
    result.put("summaries", summaries);
    result.put("totalGenes", totalGenes.size());
    result.put("calledGenes", calledGenes.size());
    result.put("hasCombo", hasCombo);
    result.put("hasMessages", hasMessages);
    result.put("hasMissingVariants", hasMissingVariants);
    result.put("hasUnphasedNote", hasUnphasedNote);
    result.put("summaryMessages", reportContext.getMessages().stream()
        .filter(MessageAnnotation.isMessage)
        .map(MessageAnnotation::getMessage)
        .toList());
    // Section III
    result.put("geneReports", geneReports);

    // Section II: Prescribing Recommendations
    SortedMap<String, Map<DataSource, DrugReport>> drugReports = new TreeMap<>();
    SortedMap<String, Recommendation> recommendationMap = new TreeMap<>();

    for (DataSource source : reportContext.getDrugReports().keySet()) {
      if (!m_sources.contains(source)) {
        continue;
      }
      for (DrugReport drugReport : reportContext.getDrugReports().get(source).values()) {
        // don't use drugReport.isMatch() directly because it escapes warfarin
        if (m_compact && drugReport.getGuidelines().stream().allMatch(GuidelineReport::isUncallable)) {
          continue;
        }

        Recommendation rec = recommendationMap.computeIfAbsent(drugReport.getName(), Recommendation::new);
        rec.addReport(source, drugReport);

        drugReports.computeIfAbsent(drugReport.getName(), (n) -> new HashMap<>())
            .put(source, drugReport);
      }
    }

    result.put("recommendations", new TreeSet<>(recommendationMap.values()));
    result.put("drugs", recommendationMap.keySet());

    return result;
  }


  private Optional<Map<String, Object>> buildGenotypeSummary(String symbol, SortedSet<GeneReport> geneReports) {

    List<GeneReport> calledReports = geneReports.stream()
        .filter(GeneReport::isReportable)
        .toList();
    if (calledReports.size() == 0) {
      return Optional.empty();
    }

    Map<String, Object> summary = new HashMap<>();
    summary.put("symbol", symbol);

    Set<String> relatedDrugs = new TreeSet<>();
    boolean hasMessages = false;
    // CPIC gets sorted first, this will pick CPIC over DPWG
    for (GeneReport report : calledReports) {
      report.getRelatedDrugs().stream()
          .map(DrugLink::getName)
          .forEach(relatedDrugs::add);
      hasMessages = hasMessages || hasMessages(report);

      if (!summary.containsKey("diplotypes")) {
        summary.put("source", report.getPhenotypeSource());
        if (report.getCallSource() == CallSource.MATCHER) {
          if (isDpyd(symbol) && report.getMatcherComponentDiplotypes().size() > 0) {
            summary.put("showComponents", true);
            summary.put("diplotypes", report.getSourceDiplotypes().get(0));
            summary.put("componentDiplotypes", report.getMatcherComponentDiplotypes());
          } else {
            summary.put("diplotypes", report.getSourceDiplotypes());
          }
        } else {
          summary.put("diplotypes", report.getSourceDiplotypes());
        }
        summary.put("hasMissingVariants", report.isMissingVariants());
        summary.put("showUnphasedNote", showUnphasedNote(report));
        summary.put("geneReport", report);

      } else {
        // TODO(markwoon): check if functionality is different?
      }
    }

    summary.put("relatedDrugs", relatedDrugs);
    summary.put("hasMessages", hasMessages);
    return Optional.of(summary);
  }

  private static boolean showUnphasedNote(GeneReport geneReport) {
    return !geneReport.isOutsideCall() && !geneReport.isPhased() && !isDpyd(geneReport.getGeneDisplay());
  }

  private static boolean hasMessages(GeneReport geneReport) {
    return geneReport.getMessages().stream().anyMatch(MessageAnnotation.isMessage);
  }


  private static class Recommendation implements Comparable<Recommendation> {
    private final String m_drug;
    private Map<String, Object> m_cpicReport;
    private Map<String, Object> m_dpwgReport;
    private boolean m_isMultiMatch;
    private boolean m_hasInferred;
    private boolean m_hasDpydInferred;
    private final Set<String> m_messages = new LinkedHashSet<>();
    private final Set<String> m_footnotes = new LinkedHashSet<>();
    private final List<Publication> m_citations = new ArrayList<>();

    public Recommendation(String drug) {
      m_drug = drug;
    }

    void addReport(DataSource source, DrugReport report) {
      Preconditions.checkArgument(m_drug.equals(report.getName()));

      switch (source) {
        case CPIC -> {
          if (m_cpicReport != null) {
            throw new IllegalStateException("Multiple drug reports for " + report.getName() + " from " +
                source.getDisplayName());
          }
          m_cpicReport = buildReport(source, report);
        }
        case DPWG -> {
          if (m_dpwgReport != null) {
            throw new IllegalStateException("Multiple drug reports for " + report.getName() + " from " +
                source.getDisplayName());
          }
          m_dpwgReport = buildReport(source, report);;
        }
      }

      // messages
      report.getMessages().stream()
          .filter(MessageAnnotation.isMessage)
          .map(MessageAnnotation::getMessage)
          .forEach(m_messages::add);
      // footnotes
      report.getMessages().stream()
          .filter(MessageAnnotation.isFootnote)
          .map(MessageAnnotation::getMessage)
          .forEach(m_footnotes::add);
      // citations
      if (report.getCitations() != null && report.getCitations().size() > 0) {
        m_citations.addAll(report.getCitations());
      }

      for (GuidelineReport guideline : report.getGuidelines()) {
        for (AnnotationReport annotation : guideline.getAnnotations()) {
          if (annotation.getGenotypes().size() > 1) {
            m_isMultiMatch = true;
          }
          for (Genotype genotype : annotation.getGenotypes()) {
            if (genotype.isInferred()) {
              if (genotype.getDiplotypes().stream()
                  .map(Diplotype::getGene)
                  .anyMatch(g -> g.equals("DPYD"))) {
                m_hasDpydInferred = true;
              } else {
                m_hasInferred = true;
              }
              break;
            }
          }
        }
      }
    }

    private Map<String, Object> buildReport(DataSource source, DrugReport drugReport) {
      Map<String, Object> report = new LinkedHashMap<>();
      report.put("id", drugReport.getId());
      report.put("matched", drugReport.isMatched());
      report.put("source", source);
      report.put("urls", drugReport.getUrls());

      if (drugReport.isMatched()) {
        report.put("guidelines", drugReport.getGuidelines());
      } else {
        boolean uncallable = drugReport.getGuidelines().stream().allMatch(GuidelineReport::isUncallable);
        report.put("uncallable", uncallable);
        if (uncallable) {
          report.put("uncalledGenes", String.join(", ", drugReport.getGuidelines().stream()
              .flatMap(gr -> gr.getUncalledGenes().stream())
              .collect(Collectors.toCollection(TreeSet::new))));
        } else {
          report.put("unmatchedCalls", String.join(" and ", drugReport.getGuidelines().stream()
              .flatMap((gr) -> gr.getRelatedGeneReports().stream())
              .flatMap((gr) -> {
                String geneLink = "<a href=\"#" + gr.getGeneDisplay() + "\">" + gr.getGeneDisplay() + "</a>";
                return gr.getRecommendationDiplotypes().stream()
                    .map((d) -> geneLink + " " + d.printBare());
              })
              .collect(Collectors.toCollection(TreeSet::new))));
        }
      }

      return report;
    }


    public String getDrug() {
      return m_drug;
    }

    public List<Map<String, Object>> getReports() {
      List<Map<String, Object>> reports = new ArrayList<>();
      if (m_cpicReport != null) {
        reports.add(m_cpicReport);
      }
      if (m_dpwgReport != null) {
        reports.add(m_dpwgReport);
      }
      return reports;
    }

    public Map<String, Object> getCpicReport() {
      return m_cpicReport;
    }

    public Map<String, Object> getDpwgReport() {
      return m_dpwgReport;
    }

    public boolean isMultiMatch() {
      return m_isMultiMatch;
    }

    public boolean isHasInferred() {
      return m_hasInferred;
    }

    public boolean isHasDpydInferred() {
      return m_hasDpydInferred;
    }

    public Set<String> getMessages() {
      return m_messages;
    }

    public Set<String> getFootnotes() {
      return m_footnotes;
    }

    public List<Publication> getCitations() {
      return m_citations;
    }


    @Override
    public int compareTo(Recommendation o) {
      return m_drug.compareTo(o.getDrug());
    }
  }
}

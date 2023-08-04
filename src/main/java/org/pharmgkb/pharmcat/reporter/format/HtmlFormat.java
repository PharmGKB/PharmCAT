package org.pharmgkb.pharmcat.reporter.format;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;
import com.github.jknack.handlebars.Handlebars;
import com.github.jknack.handlebars.helper.StringHelpers;
import com.github.jknack.handlebars.io.ClassPathTemplateLoader;
import com.google.common.collect.HashMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;
import com.google.common.collect.SortedSetMultimap;
import com.google.common.collect.TreeMultimap;
import org.apache.commons.lang3.StringUtils;
import org.pharmgkb.pharmcat.Env;
import org.pharmgkb.pharmcat.phenotype.model.GenePhenotype;
import org.pharmgkb.pharmcat.reporter.ReportContext;
import org.pharmgkb.pharmcat.reporter.format.html.Recommendation;
import org.pharmgkb.pharmcat.reporter.handlebars.ReportHelpers;
import org.pharmgkb.pharmcat.reporter.model.DataSource;
import org.pharmgkb.pharmcat.reporter.model.MessageAnnotation;
import org.pharmgkb.pharmcat.reporter.model.result.CallSource;
import org.pharmgkb.pharmcat.reporter.model.result.DrugLink;
import org.pharmgkb.pharmcat.reporter.model.result.DrugReport;
import org.pharmgkb.pharmcat.reporter.model.result.GeneReport;
import org.pharmgkb.pharmcat.reporter.model.result.GuidelineReport;

import static org.pharmgkb.pharmcat.reporter.caller.DpydCaller.isDpyd;
import static org.pharmgkb.pharmcat.reporter.caller.Slco1b1CustomCaller.isSlco1b1;


/**
 * An HTML-formatted version of {@link ReportContext} data.
 */
public class HtmlFormat extends AbstractFormat {
  private static final String sf_templatePrefix = "/org/pharmgkb/pharmcat/reporter";
  private static final String sf_handlebarTemplateName = "report";
  private static final boolean sf_compact_drugs = false;
  private List<DataSource> m_sources = Lists.newArrayList(DataSource.CPIC, DataSource.DPWG);
  private boolean m_compact;
  private final boolean f_testMode;

  /**
   * Constructor. Takes the path to write to and whether this output is "test mode" or not
   * @param outputPath the path to write the data to
   * @param testMode if true will leave out some time-related metadata that may give false-positives while diffing
   * output
   */
  public HtmlFormat(Path outputPath, Env env, boolean testMode) {
    super(outputPath, env);
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
      writer.write(handlebars.compile(sf_handlebarTemplateName).apply(reportData));
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
      result.put("timestamp", new SimpleDateFormat("MMMM dd, yyyy").format(reportContext.timestamp()));
      result.put("pharmcatVersion", reportContext.getPharmcatVersion());
      result.put("cpicVersion", reportContext.getCpicVersion());
    }
    result.put("compact", m_compact);

    if (StringUtils.isNotBlank(reportContext.getTitle())) {
      result.put("title", reportContext.getTitle());
    }


    // Section I: Genotype Summary
    SortedSet<String> totalGenes = new TreeSet<>();
    SortedSet<String> calledGenes = new TreeSet<>();
    SortedSet<String> noDataGenes = new TreeSet<>();
    SortedSet<String> uncallableGenes = new TreeSet<>();
    boolean hasCombo = false;
    boolean hasMessages = false;
    boolean hasMissingVariants = false;
    boolean hasUnphasedNote = false;
    boolean hasUndocumentedVariationsAsReference = false;

    SortedSetMultimap<String, GeneReport> geneReportMap = TreeMultimap.create();
    Map<String, Map<String, String>> functionMap = new HashMap<>();
    for (DataSource source : new TreeSet<>(reportContext.getGeneReports().keySet())) {
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

        // CPIC gets sorted first, this will pick CPIC over DPWG
        if (!functionMap.containsKey(symbol)) {
          GenePhenotype genePhenotype = getEnv().getPhenotype(symbol, source);
          if (genePhenotype != null) {
            functionMap.put(symbol, genePhenotype.getHaplotypes());
          }
        }

        if (geneReport.isNoData()) {
          noDataGenes.add(symbol);
          continue;
        }

        // skip gene reports that aren't related to any drugs
        if (geneReport.getRelatedDrugs().size() == 0) {
          continue;
        }

        if (m_compact) {
          geneReportMap.put(symbol, geneReport);
        }
        if (geneReport.isReportable()) {
          calledGenes.add(symbol);

          if (geneReport.isMissingVariants()) {
            hasMissingVariants = true;
          }
          hasCombo = hasCombo || geneReport.getMessages().stream()
              .anyMatch(m -> m.getExceptionType().equals(MessageAnnotation.TYPE_COMBO));
          hasMessages = hasMessages || hasMessages(geneReport);
          hasUnphasedNote = hasUnphasedNote || showUnphasedNote(geneReport);
        } else {
          uncallableGenes.add(symbol);
        }
      }
    }

    SortedSet<String> genes = new TreeSet<>(geneReportMap.keySet());
    List<Map<String, Object>> summaries = new ArrayList<>();
    List<GeneReport> geneReports = new ArrayList<>();
    Multimap<String, Map<String, Object>> geneSummariesByDrug = HashMultimap.create();
    for (String symbol : genes) {
      if (noDataGenes.contains(symbol)) {
        continue;
      }
      SortedSet<GeneReport> reports = geneReportMap.get(symbol);
      if (reports.size() > 2) {
        throw new IllegalStateException("More than 2 gene reports for " + symbol);
      }
      if (reports.stream().noneMatch(GeneReport::isReportable)) {
        // add to gene reports to show in Section III
        // if (a) extended mode or (b) cannot be called but has VCF data
        if (!m_compact || uncallableGenes.contains(symbol) ||
            reports.stream().anyMatch(GeneReport::isHasUndocumentedVariations)) {
          geneReports.add(reports.first());
        }
        continue;
      }
      if (reports.first().isTreatUndocumentedVariationsAsReference()) {
        hasUndocumentedVariationsAsReference = true;
      }
      Map<String, Object> geneSummary = buildGenotypeSummary(symbol, reports);
      summaries.add(geneSummary);
      //noinspection unchecked
      ((Set<String>)geneSummary.get("relatedDrugs"))
          .forEach(d -> geneSummariesByDrug.put(d, geneSummary));
      geneReports.add(reports.first());
    }
    result.put("genes", genes);
    result.put("noDataGenes", noDataGenes);
    result.put("uncallableGenes", uncallableGenes);
    result.put("summaries", summaries);
    result.put("totalGenes", totalGenes.size());
    result.put("calledGenes", calledGenes.size());
    result.put("hasCombo", hasCombo);
    result.put("hasMessages", hasMessages);
    result.put("hasMissingVariants", hasMissingVariants);
    result.put("hasUndocumentedVariationsAsReference", hasUndocumentedVariationsAsReference);
    result.put("hasUnphasedNote", hasUnphasedNote);
    result.put("summaryMessages", reportContext.getMessages().stream()
        .filter(MessageAnnotation.isMessage)
        .toList());
    // Section III
    result.put("geneReports", geneReports);
    result.put("functionMap", functionMap);

    // Section II: Prescribing Recommendations
    SortedMap<String, Map<DataSource, DrugReport>> drugReports = new TreeMap<>();
    SortedMap<String, Recommendation> recommendationMap = new TreeMap<>();

    for (DataSource source : reportContext.getDrugReports().keySet()) {
      if (!m_sources.contains(source)) {
        continue;
      }
      for (DrugReport drugReport : reportContext.getDrugReports().get(source).values()) {
        // don't use drugReport.isMatch() directly because it escapes warfarin
        if (m_compact && drugReport.getGuidelines().stream().noneMatch(GuidelineReport::isReportable)) {
          continue;
        }

        Recommendation rec = recommendationMap.computeIfAbsent(drugReport.getName(), n -> new Recommendation(getEnv(), n));
        rec.addReport(source, drugReport);

        drugReports.computeIfAbsent(drugReport.getName(), (n) -> new HashMap<>())
            .put(source, drugReport);
      }
    }

    SortedSet<Recommendation> recommendations = new TreeSet<>();
    SortedSet<String> drugsWithRecommendations = new TreeSet<>();
    if (m_compact && sf_compact_drugs) {
      for (Recommendation recommendation : recommendationMap.values()) {
        if (recommendation.isMatched()) {
          recommendations.add(recommendation);
          drugsWithRecommendations.add(recommendation.getDrug());
        }
      }
    } else {
      recommendations.addAll(recommendationMap.values());
      drugsWithRecommendations.addAll(recommendationMap.keySet());
    }

    result.put("recommendations", recommendations);
    result.put("drugsWithRecommendations", drugsWithRecommendations);
    result.put("drugs", recommendationMap.keySet());

    return result;
  }


  private Map<String, Object> buildGenotypeSummary(String symbol, SortedSet<GeneReport> geneReports) {

    Map<String, Object> summary = new HashMap<>();
    summary.put("symbol", symbol);

    Set<String> relatedDrugs = new TreeSet<>();
    boolean hasMessages = false;
    // CPIC gets sorted first, this will pick CPIC over DPWG
    for (GeneReport report : geneReports) {
      if (!report.isReportable()) {
        continue;
      }
      report.getRelatedDrugs().stream()
          .map(DrugLink::getName)
          .forEach(relatedDrugs::add);
      hasMessages = hasMessages || hasMessages(report);

      if (!summary.containsKey("diplotypes")) {
        summary.put("source", report.getPhenotypeSource());
        if (report.getCallSource() == CallSource.MATCHER) {
          if (isSlco1b1(symbol)) {
            summary.put("diplotypes", report.getRecommendationDiplotypes());
          } else if (isDpyd(symbol) && !report.getMatcherComponentHaplotypes().isEmpty()) {
            summary.put("showComponents", true);
            summary.put("diplotypes", report.getSourceDiplotypes().first());
            summary.put("componentDiplotypes", report.getMatcherComponentHaplotypes());
          } else {
            summary.put("diplotypes", report.getSourceDiplotypes());
          }
        } else {
          summary.put("diplotypes", report.getRecommendationDiplotypes());
        }
        summary.put("homozygousComponentHaplotypes", report.getMatcherHomozygousComponentHaplotypes());
        summary.put("hasMissingVariants", report.isMissingVariants());
        summary.put("showUnphasedNote", showUnphasedNote(report));
        summary.put("hasUndocumentedVariants", report.isHasUndocumentedVariations());
        summary.put("treatUndocumentedVariationsAsReference", report.isTreatUndocumentedVariationsAsReference());
        summary.put("geneReport", report);

      } else {
        // TODO(markwoon): do we need to do anything special when there are GeneReports from more than 1 source?
      }
    }

    summary.put("relatedDrugs", relatedDrugs);
    summary.put("hasMessages", hasMessages);
    return summary;
  }

  private static boolean showUnphasedNote(GeneReport geneReport) {
    return !geneReport.isOutsideCall() && !geneReport.isPhased() && !isDpyd(geneReport.getGeneDisplay());
  }

  private static boolean hasMessages(GeneReport geneReport) {
    return geneReport.getMessages().stream().anyMatch(MessageAnnotation.isMessage);
  }
}

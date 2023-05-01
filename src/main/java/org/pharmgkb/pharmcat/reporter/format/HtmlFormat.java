package org.pharmgkb.pharmcat.reporter.format;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.stream.Collectors;
import com.github.jknack.handlebars.Handlebars;
import com.github.jknack.handlebars.helper.StringHelpers;
import com.github.jknack.handlebars.io.ClassPathTemplateLoader;
import com.google.common.base.Preconditions;
import com.google.common.collect.HashMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;
import com.google.common.collect.SortedSetMultimap;
import com.google.common.collect.TreeMultimap;
import org.apache.commons.lang3.StringUtils;
import org.pharmgkb.pharmcat.Env;
import org.pharmgkb.pharmcat.phenotype.model.GenePhenotype;
import org.pharmgkb.pharmcat.reporter.MessageHelper;
import org.pharmgkb.pharmcat.reporter.ReportContext;
import org.pharmgkb.pharmcat.reporter.handlebars.ReportHelpers;
import org.pharmgkb.pharmcat.reporter.model.DataSource;
import org.pharmgkb.pharmcat.reporter.model.MessageAnnotation;
import org.pharmgkb.pharmcat.reporter.model.cpic.Publication;
import org.pharmgkb.pharmcat.reporter.model.result.AnnotationReport;
import org.pharmgkb.pharmcat.reporter.model.result.CallSource;
import org.pharmgkb.pharmcat.reporter.model.result.Diplotype;
import org.pharmgkb.pharmcat.reporter.model.result.DrugLink;
import org.pharmgkb.pharmcat.reporter.model.result.DrugReport;
import org.pharmgkb.pharmcat.reporter.model.result.GeneReport;
import org.pharmgkb.pharmcat.reporter.model.result.Genotype;
import org.pharmgkb.pharmcat.reporter.model.result.GuidelineReport;

import static org.pharmgkb.pharmcat.reporter.caller.DpydCaller.isDpyd;


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
        if (!m_compact) {
          geneReports.add(reports.first());
        }
        continue;
      }
      if (reports.first().isTreatUndocumentedVariationsAsReference()) {
        hasUndocumentedVariationsAsReference = true;
      }
      summaries.add(buildGenotypeSummary(symbol, reports));
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
          if (isDpyd(symbol) && report.getMatcherComponentHaplotypes().size() > 0) {
            summary.put("showComponents", true);
            summary.put("diplotypes", report.getSourceDiplotypes().get(0));
            summary.put("componentDiplotypes", report.getMatcherComponentHaplotypes());
          } else {
            summary.put("diplotypes", report.getSourceDiplotypes());
          }
        } else {
          summary.put("diplotypes", report.getSourceDiplotypes());
        }
        summary.put("hasMissingVariants", report.isMissingVariants());
        summary.put("showUnphasedNote", showUnphasedNote(report));
        summary.put("hasUndocumentedVariants", report.isHasUndocumentedVariations());
        summary.put("treatUndocumentedVariationsAsReference", report.isTreatUndocumentedVariationsAsReference());
        summary.put("geneReport", report);

      } else {
        // TODO(markwoon): check if functionality is different?
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


  private class Recommendation implements Comparable<Recommendation> {
    private final String m_drug;
    private Map<String, Object> m_cpicReport;
    private Map<String, Object> m_dpwgReport;
    private boolean m_hasInferred;
    private boolean m_hasDpydInferred;
    private final Set<MessageAnnotation> m_messages = new LinkedHashSet<>();
    private final Set<MessageAnnotation> m_footnotes = new LinkedHashSet<>();
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
          m_dpwgReport = buildReport(source, report);
        }
      }

      boolean callMultiMatch = false;
      boolean scoreMultimatch = false;
      Multimap<String, Diplotype> mismatchDiplotypes = HashMultimap.create();
      for (GuidelineReport guideline : report.getGuidelines()) {
        for (AnnotationReport annotation : guideline.getAnnotations()) {
          if (annotation.getGenotypes().size() > 1) {
            callMultiMatch = true;
          }
          for (Genotype genotype : annotation.getGenotypes()) {
            for (Diplotype diplotype : genotype.getDiplotypes()) {
              if (!scoreMultimatch && diplotype.getLookupKeys().size() > 1) {
                GenePhenotype gp = getEnv().getPhenotype(diplotype.getGene(), source);
                scoreMultimatch = gp != null && gp.isMatchedByActivityScore();
              }
              if (diplotype.getOutsidePhenotypeMismatch() != null ||
                  diplotype.getOutsideActivityScoreMismatch() != null) {
                mismatchDiplotypes.put(diplotype.getGene(), diplotype);
              }
            }
            if (scoreMultimatch || genotype.getDiplotypes().stream().anyMatch((d) -> {
              if (d.getLookupKeys().size() > 1) {
                GenePhenotype gp = getEnv().getPhenotype(d.getGene(), source);
                return gp != null && gp.isMatchedByActivityScore();
              }
              return false;
            })) {
              scoreMultimatch = true;
            }
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
      if (callMultiMatch) {
        m_messages.add(getEnv().getMessage(MessageHelper.MSG_MUlTI_CALL));
      }
      if (scoreMultimatch) {
        m_messages.add(getEnv().getMessage(MessageHelper.MSG_MULTI_SCORE));
      }
      for (String gene : mismatchDiplotypes.keySet()) {
        StringBuilder builder = new StringBuilder()
            .append("Conflicting outside call data was provided for ")
            .append(gene)
            .append(".  PharmCAT will use provided ");
        if (mismatchDiplotypes.get(gene).stream().anyMatch(Diplotype::hasActivityScore)) {
          builder.append("activity score");
        } else {
          builder.append("phenotype");
        }
        builder.append(" to match recommendations.");
        m_messages.add(new MessageAnnotation(MessageAnnotation.TYPE_NOTE, "warn.mismatch.outsideCalll",
            builder.toString()));
      }

      // messages
      report.getMessages().stream()
          .filter(MessageAnnotation.isMessage)
          .forEach(m_messages::add);
      // footnotes
      report.getMessages().stream()
          .filter(MessageAnnotation.isFootnote)
          .forEach(m_footnotes::add);
      // citations
      if (report.getCitations() != null && report.getCitations().size() > 0) {
        m_citations.addAll(report.getCitations());
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
              .flatMap((guidelineReport) -> guidelineReport.getGeneReports().stream())
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

    public boolean isHasInferred() {
      return m_hasInferred;
    }

    public boolean isHasDpydInferred() {
      return m_hasDpydInferred;
    }

    public Set<MessageAnnotation> getMessages() {
      return m_messages;
    }

    public Set<MessageAnnotation> getFootnotes() {
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

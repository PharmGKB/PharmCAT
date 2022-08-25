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
import java.util.List;
import java.util.Map;
import java.util.TreeSet;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import com.github.jknack.handlebars.Handlebars;
import com.github.jknack.handlebars.helper.StringHelpers;
import com.github.jknack.handlebars.io.ClassPathTemplateLoader;
import org.apache.commons.lang3.StringUtils;
import org.pharmgkb.pharmcat.reporter.ReportContext;
import org.pharmgkb.pharmcat.reporter.TextConstants;
import org.pharmgkb.pharmcat.reporter.handlebars.ReportHelpers;
import org.pharmgkb.pharmcat.reporter.model.MessageAnnotation;
import org.pharmgkb.pharmcat.reporter.model.VariantReport;
import org.pharmgkb.pharmcat.reporter.model.result.Diplotype;
import org.pharmgkb.pharmcat.reporter.model.result.DrugReport;
import org.pharmgkb.pharmcat.reporter.model.result.GeneReport;

import static org.pharmgkb.pharmcat.reporter.TextConstants.NA;
import static org.pharmgkb.pharmcat.reporter.model.result.GeneReport.UNCALLED;


/**
 * An HTML-formatted version of {@link ReportContext} data.
 */
public class HtmlFormat extends AbstractFormat {
  private static final String sf_templatePrefix = "/org/pharmgkb/pharmcat/reporter";
  private static final String FINAL_REPORT = "report";
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
   * handlebars template.
   * @return a Map of data to serialize into JSON
   */
  private Map<String,Object> compile(ReportContext reportContext) {

    Map<String,Object> result = new HashMap<>();
    if (!f_testMode) {
      result.put("generatedOn", new SimpleDateFormat("MMMMM dd, yyyy").format(reportContext.getGeneratedOn()));
      result.put("pharmcatVersion", reportContext.getPharmcatVersion());
    }

    if (StringUtils.isNotBlank(reportContext.getTitle())) {
      result.put("title", reportContext.getTitle());
    }

    // Genotypes section
    List<Map<String,Object>> genotypes = new ArrayList<>();
    int calledGenes = 0;
    int totalGenes = 0;
    boolean hasCombo = false;
    for (GeneReport geneReport : reportContext.getGeneReports()) {
      // skip any genes on the blacklist
      if (geneReport.isIgnored()) {
        continue;
      }
      totalGenes += 1;

      // skip any uncalled genes
      boolean allVariantsMissing = geneReport.getVariantReports().stream().allMatch(VariantReport::isMissing);
      if ((!geneReport.isCalled() || allVariantsMissing) && (geneReport.getReporterDiplotypes().isEmpty())) {
        continue;
      }

      if (geneReport.getRelatedDrugs().size() == 0) {
        continue;
      }

      Map<String,Object> genotype = new HashMap<>();
      genotype.put("gene", geneReport.getGeneDisplay());
      genotype.put("called", geneReport.isCalled());
      genotype.put("reportable", geneReport.isReportable());
      genotype.put("drugs", geneReport.getRelatedDrugs());
      genotype.put("diplotypes", makeSummaryDiplotypes(geneReport));
      genotype.put("missingVariants", geneReport.isMissingVariants());
      genotype.put("unphased", !geneReport.isOutsideCall() && !geneReport.isPhased());
      genotype.put("hasMessages", geneReport.getMessages().size()>0);
      genotype.put("outsideCall", geneReport.isOutsideCall());

      genotypes.add(genotype);

      hasCombo = hasCombo || geneReport.getMessages().stream().anyMatch(m -> m.getExceptionType().equals(MessageAnnotation.TYPE_COMBO));

      if (geneReport.isReportable()) {
        calledGenes += 1;
      }
    }
    result.put("genotypes", genotypes);
    result.put("totalGenes", totalGenes);
    result.put("calledGenes", calledGenes);
    result.put("hasCombo", hasCombo);
    result.put("messages", reportContext.getMessages().stream().map(MessageAnnotation::getMessage).collect(Collectors.toList()));

    // Drugs section
    List<Map<String,Object>> drugReports = new ArrayList<>();
    for (DrugReport drugReport : reportContext.getDrugReports()) {
      Map<String,Object> drugMap = new LinkedHashMap<>();

      drugMap.put("id", drugReport.getId());
      drugMap.put("name", drugReport.getName());
      drugMap.put("urls", drugReport.getUrls());

      List<Map<String,Object>> geneCallList = new ArrayList<>();
      for (String gene : drugReport.getRelatedGeneSymbols()) {
        reportContext.findGeneReport(gene).ifPresent((geneReport) -> {
          String functions = geneReport.isReportable() ? String.join("; ", geneReport.printDisplayFunctions()) : null;
          Map<String,Object> geneCall = new LinkedHashMap<>();
          geneCall.put("gene", geneReport.getGeneDisplay());
          geneCall.put("diplotypes", String.join(", ", geneReport.printDisplayCalls()));
          geneCall.put("showHighlights", !geneReport.getHighlightedVariants().isEmpty());
          geneCall.put("highlightedVariants", geneReport.getHighlightedVariants());
          geneCall.put("functions", functions);
          geneCall.put("outsideCall", geneReport.isOutsideCall());
          geneCallList.add(geneCall);
        });
      }
      for (String variant : drugReport.getReportVariants()) {
        String call = reportContext.getGeneReports().stream()
            .flatMap(g -> Stream.concat(g.getVariantReports().stream(), g.getVariantOfInterestReports().stream()))
            .filter(v -> v.getDbSnpId() != null && v.getDbSnpId().matches(variant) && !v.isMissing())
            .map(VariantReport::getCall)
            .collect(Collectors.joining(", "));
        Map<String,Object> geneCall = new LinkedHashMap<>();
        geneCall.put("gene", variant);
        geneCall.put("diplotypes", StringUtils.isBlank(call) ? "missing" : call);
        geneCall.put("outsideCall", false);
        geneCallList.add(geneCall);
      }
      if (geneCallList.size() > 0) {
        drugMap.put("geneCalls", geneCallList);
      }

      drugMap.put("matched", drugReport.isMatched());

      drugMap.put("messages", drugReport.getMessages().stream()
          .filter(MessageAnnotation.isMessage)
          .map(MessageAnnotation::getMessage)
          .collect(Collectors.toList()));

      drugMap.put("footnotes", drugReport.getMessages().stream()
          .filter(MessageAnnotation.isFootnote)
          .map(MessageAnnotation::getMessage)
          .collect(Collectors.toList()));

      if (drugReport.getCitations() != null && drugReport.getCitations().size()>0) {
        drugMap.put("citations", drugReport.getCitations());
      }

      // special case the display for warfarin recommendation since it's an image
      if (drugReport.toString().equals("warfarin")) {
        Map<String,String> imageData = new LinkedHashMap<>();
        imageData.put("url", "https://files.cpicpgx.org/images/warfarin/warfarin_recommendation_diagram.png");
        imageData.put("altText", "Figure 2 from the CPIC guideline for warfarin");
        drugMap.put("image", imageData);
      }

      drugMap.put("guidelines", drugReport.getGuidelines());

      drugReports.add(drugMap);
    }
    result.put("guidelines", drugReports);


    // Gene calls

    List<Map<String,Object>> geneCallList = new ArrayList<>();
    for (GeneReport geneReport : reportContext.getGeneReports()) {
      if (geneReport.isIgnored()) {
        continue;
      }

      Map<String,Object> geneCallMap = new HashMap<>();

      geneCallMap.put("gene", geneReport.getGeneDisplay());

      String phaseStatus;
      boolean unphased = false;
      if (geneReport.isOutsideCall()) {
        phaseStatus = "Unavailable for calls made outside PharmCAT";
      } else {
        phaseStatus = geneReport.isPhased() ? "Phased" : "Unphased";
        unphased = !geneReport.isPhased();
      }
      geneCallMap.put("phaseStatus", phaseStatus);
      geneCallMap.put("unphased", unphased);

      boolean hasUncalledHaplotypes = geneReport.getUncalledHaplotypes() != null && geneReport.getUncalledHaplotypes().size() > 0;
      geneCallMap.put("hasUncalledHaps", hasUncalledHaplotypes);
      if (hasUncalledHaplotypes) {
        geneCallMap.put("uncalledHaps", String.join(", ", geneReport.getUncalledHaplotypes()));
      }

      List<String> diplotypes = new ArrayList<>(geneReport.printDisplayCalls());
      diplotypes.addAll(geneReport.getHighlightedVariants());
      geneCallMap.put("diplotypes", diplotypes);

      if (geneReport.getMessages() != null && geneReport.getMessages().size() > 0) {
        geneCallMap.put("warnings", geneReport.getMessages().stream().map(MessageAnnotation::getMessage).collect(Collectors.toList()));
      }

      if (geneReport.getVariantReports().size() > 0) {
        geneCallMap.put("variants", new TreeSet<>(geneReport.getVariantReports()));
      } else {
        geneCallMap.put("variantsUnspecified", true);
      }

      if (geneReport.getVariantOfInterestReports().size() > 0) {
        geneCallMap.put("variantsOfInterest", geneReport.getVariantOfInterestReports());
      } else {
        geneCallMap.put("variantsOfInterestUnspecified", true);
      }

      geneCallMap.put("outsideCall", geneReport.isOutsideCall());

      geneCallMap.put("totalMissingVariants",
          geneReport.getVariantReports().stream().filter(VariantReport::isMissing).count());
      geneCallMap.put("totalVariants", geneReport.getVariantReports().size());

      geneCallMap.put("messages", geneReport.getMessages().stream()
          .filter(MessageAnnotation.isMessage)
          .map(MessageAnnotation::getMessage)
          .collect(Collectors.toList()));
      geneCallMap.put("extra-position-notes", geneReport.getMessages().stream()
          .filter(MessageAnnotation.isExtraPositionNote)
          .map(MessageAnnotation::getMessage)
          .collect(Collectors.toList()));

      geneCallList.add(geneCallMap);
    }
    result.put("geneCalls", geneCallList);

    return result;
  }

  private List<HtmlDiplotype> makeSummaryDiplotypes(GeneReport geneReport) {
    List<Diplotype> diplotypesToUse = geneReport.isLeastFunction()
        ? geneReport.getMatcherDiplotypes()
        : geneReport.getReporterDiplotypes();
    return diplotypesToUse.stream().map(HtmlDiplotype::new).toList();
  }

  private static class HtmlDiplotype {
    private String m_call = UNCALLED;
    private String m_function;
    private String m_phenotype = NA;

    HtmlDiplotype(Diplotype diplotype) {
      if (diplotype == null) {
        return;
      }
      m_call = diplotype.printDisplay();
      if (diplotype.isCombination() && diplotype.getGene().equals("DPYD")) {
        m_function = TextConstants.SEE_DRUG;
      } else {
        m_function = diplotype.printFunctionPhrase();
      }
      m_phenotype = diplotype.printPhenotypes();
    }

    public String getCall() {
      return m_call;
    }

    public String getFunction() {
      return m_function;
    }

    public String getPhenotype() {
      return m_phenotype;
    }
  }
}

package org.pharmgkb.pharmcat.reporter;

import java.io.IOException;
import java.lang.invoke.MethodHandles;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Date;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import javax.annotation.Nonnull;
import org.apache.commons.lang3.StringUtils;
import org.checkerframework.checker.nullness.qual.Nullable;
import org.pharmgkb.pharmcat.ParseException;
import org.pharmgkb.pharmcat.definition.PhenotypeMap;
import org.pharmgkb.pharmcat.haplotype.DefinitionReader;
import org.pharmgkb.pharmcat.haplotype.model.GeneCall;
import org.pharmgkb.pharmcat.reporter.model.MessageAnnotation;
import org.pharmgkb.pharmcat.reporter.model.OutsideCall;
import org.pharmgkb.pharmcat.reporter.model.VariantReport;
import org.pharmgkb.pharmcat.reporter.model.cpic.Drug;
import org.pharmgkb.pharmcat.reporter.model.cpic.Recommendation;
import org.pharmgkb.pharmcat.reporter.model.result.CallSource;
import org.pharmgkb.pharmcat.reporter.model.result.DrugReport;
import org.pharmgkb.pharmcat.reporter.model.result.GeneReport;
import org.pharmgkb.pharmcat.util.CliUtils;
import org.pharmgkb.pharmcat.util.DataManager;
import org.pharmgkb.pharmcat.util.MessageMatcher;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


/**
 * This class acts as a central context for all data needed to generate the final report.
 *
 * It currently gathers
 * <ul>
 *   <li>{@link GeneCall} objects from the named allele matcher</li>
 *   <li>{@link DrugReport} objects from dosing guideline annotations</li>
 *   <li>Allele definitions on a per-gene basis</li>
 * </ul>
 *
 * @author greytwist
 * @author Ryan Whaley
 */
public class ReportContext {
  private static final Logger sf_logger = LoggerFactory.getLogger(MethodHandles.lookup().lookupClass());

  private final Map<String,GeneReport> m_geneReports = new TreeMap<>();
  private final List<DrugReport> m_drugReports = new ArrayList<>();
  private final PhenotypeMap m_phenotypeMap;
  private final Map<String,String> m_refAlleleForGene = new HashMap<>();

  /**
   * Public constructor. Compiles all the incoming data into useful objects to be held for later reporting
   * @param calls {@link GeneCall} objects from the sample data
   * @param outsideCalls {@link OutsideCall} objects, non-null but can be empty
   */
  public ReportContext(List<GeneCall> calls, List<OutsideCall> outsideCalls) throws Exception {
    for (Drug drug : new DrugCollection()) {
      DrugReport drugReport = new DrugReport(drug);

      // do not report any drug (or gene) that has a gene on the "ignored" list
      if (drugReport.isIgnored()) continue;

      // initialize the gene report
      m_drugReports.add(drugReport);

      // initialize any non-present drug reports
      for (String gene : drug.getGenes()) {
        m_geneReports.putIfAbsent(gene, new GeneReport(gene));
      }
    }

    loadReferenceAlleleNames();

    m_phenotypeMap = new PhenotypeMap();

    compileMatcherData(calls);
    compileOutsideCallData(outsideCalls);

    findMatches();

    compileDrugsForGenes();
  }

  /**
   * Applies the given {@link MessageAnnotation} objects to the data in this report.
   * @param messages a List of {@link MessageAnnotation} objects
   */
  public void applyMessage(List<MessageAnnotation> messages) {
    MessageMatcher messageMatcher = new MessageMatcher(messages, this);

    m_geneReports.values().forEach(r -> r.addMessages(messageMatcher.match(r)));
    m_drugReports.forEach(r -> r.addMessages(messageMatcher.match(r)));
  }

  /**
   * Takes {@link GeneCall} data and preps internal data structures for usage. Also prepares exception logic and applies
   * it to the calling data
   */
  private void compileMatcherData(List<GeneCall> calls) throws Exception {
    for (GeneCall call : calls) {
      GeneReport geneReport = new GeneReport(call.getGene());
      geneReport.setCallData(call);
      m_geneReports.put(call.getGene(), geneReport);

      DiplotypeFactory diplotypeFactory = new DiplotypeFactory(
          call.getGene(),
          m_phenotypeMap.lookup(call.getGene()).orElse(null),
          m_refAlleleForGene.get(call.getGene()));
      geneReport.setDiplotypes(diplotypeFactory, call);
    }
  }

  /**
   * Takes outside calls, find the GeneReport for each one and then adds call information to it
   * @param calls outside calls
   */
  private void compileOutsideCallData(List<OutsideCall> calls) {
    for (OutsideCall outsideCall : calls) {
      GeneReport geneReport = getGeneReport(outsideCall.getGene());
      if (geneReport.isOutsideCall()) {
        throw new ParseException("Duplicate outside call found for " + geneReport.getGene());
      }
      if (geneReport.getCallSource() == CallSource.MATCHER && geneReport.isCalled()) {
        throw new ParseException("Cannot specify outside call for " + geneReport.getGene() + " since it's already in sample data");
      }

      geneReport.setOutsideCallData();

      DiplotypeFactory diplotypeFactory = new DiplotypeFactory(
          outsideCall.getGene(),
          m_phenotypeMap.lookup(outsideCall.getGene()).orElse(null),
          m_refAlleleForGene.get(outsideCall.getGene()));
      geneReport.setDiplotypes(diplotypeFactory, outsideCall);
    }
  }

  /**
   * Adds compiled drug data to each {@link GeneReport} object that's related through a guideline. This gives a useful
   * collection of all drugs that are related to each gene.
   */
  private void compileDrugsForGenes() {
    m_drugReports.forEach(r -> {
      for (String gene : r.getRelatedGeneSymbols()) {
        getGeneReport(gene).addRelatedDrugs(r);
      }
    });
  }

  private boolean isReportable(String gene) {
    return getGeneReport(gene).isReportable();
  }

  /**
   * Assigns matched guideline groups for all guidelines in this report based on called diplotype functions for each
   * gene and guideline combination.
   */
  private void findMatches() {

    for(DrugReport drugReport : m_drugReports) {
      boolean reportable = drugReport.getRelatedGeneSymbols().stream()
          .anyMatch(this::isReportable);

      drugReport.setReportable(reportable);

      drugReport.getRelatedGeneSymbols().stream()
          .filter(g -> !isReportable(g))
          .forEach(drugReport::addUncalledGene);

      if (!reportable) {
        continue;
      }

      makeAllReportGenotypes(drugReport);
    }
  }

  /**
   *
   */
  private void loadReferenceAlleleNames() throws IOException {
    m_refAlleleForGene.put("CYP2D6", "*1");

    DefinitionReader definitionReader = new DefinitionReader();
    definitionReader.read(DataManager.DEFAULT_DEFINITION_DIR);
    definitionReader.getGenes()
        .forEach(g -> m_refAlleleForGene.put(g, definitionReader.getHaplotypes(g).get(0).getName()));
  }

  /**
   * Makes a set of called genotype function Strings (in the form "GENEA:No Function/No Function;GENEB:Normal Function/Normal Function") for the given
   * {@link DrugReport} and then adds them to the given {@link DrugReport}.
   *
   * <em>note:</em> this needs to stay in the {@link ReportContext} since it relies on referencing possibly multiple
   * {@link GeneReport} objects.
   */
  private void makeAllReportGenotypes(DrugReport drugReport) {
    Set<String> results = new TreeSet<>();
    for (String symbol : drugReport.getRelatedGeneSymbols()) {
      results = makeCalledGenotypes(symbol, results);
    }
    results.forEach(drugReport::addReportGenotype);
  }

  /**
   * Makes the genotype Strings for the given gene <code>symbol</code> and either creates a new Set if the passed
   * <code>results</code> is empty or adds to the existing <code>results</code> if there are already entries
   *
   * @param symbol the gene to generate calls for
   * @param results the existing call list to add to
   * @return a new Set of gene calls with the calls for the specified gene
   */
  private Set<String> makeCalledGenotypes(String symbol, Set<String> results) {
    if (results.size() == 0) {
      return getGeneReport(symbol).getDiplotypeLookupKeys().stream()
          .map(k -> m_phenotypeMap.lookupPhenotype(k).orElse(symbol+":N/A"))
          .collect(Collectors.toSet());
    }
    else {
      Set<String> newResults = new TreeSet<>();
      for (String geno1 : results) {
        getGeneReport(symbol).getDiplotypeLookupKeys().stream().map(m_phenotypeMap::lookupPhenotype).forEach(
            geno2 -> {
              Set<String> genos = new TreeSet<>();
              genos.add(geno1);
              genos.add(geno2.orElse(symbol+":N/A"));
              newResults.add(String.join(";", genos));
            });
      }
      return newResults;
    }
  }

  public List<DrugReport> getDrugReports() {
    return m_drugReports;
  }

  public Collection<GeneReport> getGeneReports() {
    return m_geneReports.values();
  }

  @Nonnull
  public GeneReport getGeneReport(String geneSymbol) {
    return m_geneReports.get(geneSymbol);
  }

  /**
   * Make a Map of data that will be used in the final report. This map will be serialized and then applied to the 
   * handlebars template.
   * @return a Map of data to serialize into JSON
   */
  public Map<String,Object> compile(@Nullable String title) throws IOException {

    Map<String,Object> result = new HashMap<>();
    result.put("generatedOn", new SimpleDateFormat("MMMMM dd, yyyy").format(new Date()));
    result.put("version", CliUtils.getVersion());

    if (StringUtils.isNotBlank(title)) {
      result.put("title", title);
    }

    // Genotypes section
    List<Map<String,Object>> genotypes = new ArrayList<>();
    int calledGenes = 0;
    for (GeneReport geneReport : getGeneReports()) {
      String symbol = geneReport.getGene();

      // skip any genes on the blacklist
      if (geneReport.isIgnored()) {
        continue;
      }

      // skip any uncalled genes
      if (!geneReport.isCalled() && geneReport.getReporterDiplotypes().isEmpty()) {
        continue;
      }

      if (geneReport.getRelatedDrugs().size() == 0) {
        continue;
      }

      Map<String,Object> genotype = new HashMap<>();
      genotype.put("gene", symbol);
      genotype.put("called", geneReport.isCalled());
      genotype.put("drugs", geneReport.getRelatedDrugs());
      genotype.put("calls", geneReport.printDisplayCalls());
      genotype.put("functions", geneReport.printDisplayFunctions());
      genotype.put("missingVariants", geneReport.isMissingVariants());
      genotype.put("phenotype", geneReport.printDisplayPhenotypes());
      genotype.put("hasMessages", geneReport.getMessages().size()>0);
      genotype.put("outsideCall", geneReport.isOutsideCall());

      genotypes.add(genotype);

      if (geneReport.isCalled()) {
        calledGenes += 1;
      }
    }
    result.put("genotypes", genotypes);
    result.put("totalGenes", genotypes.size());
    result.put("calledGenes", calledGenes);

    // Guidelines section
    List<Map<String,Object>> guidelines = new ArrayList<>();
    for (DrugReport guideline : getDrugReports()) {
      sf_logger.debug("Start GuidelineReport for {}", guideline.toString());

      Map<String,Object> guidelineMap = new HashMap<>();

      String drugs = String.join(", ", guideline.getRelatedDrugs());
      guidelineMap.put("drugs", drugs);
      guidelineMap.put("url", guideline.getUrl());
      guidelineMap.put("id", guideline.getId());
      String lastModified = guideline.getLastModified() != null
          ? new SimpleDateFormat("yyyy.MM.dd").format(guideline.getLastModified())
          : "not available";
      guidelineMap.put("lastModified", lastModified);

      List<Map<String,Object>> geneCallList = new ArrayList<>();
      for (String gene : guideline.getRelatedGeneSymbols()) {
        GeneReport geneReport = getGeneReport(gene);
        String functions = geneReport.isCalled() ? String.join("; ", geneReport.printDisplayFunctions()) : null;
        Map<String,Object> geneCall = new LinkedHashMap<>();
        geneCall.put("gene", gene);
        geneCall.put("diplotypes", String.join(", ", geneReport.printDisplayCalls()));
        geneCall.put("showHighlights", !geneReport.getHighlightedVariants().isEmpty());
        geneCall.put("highlightedVariants", geneReport.getHighlightedVariants());
        geneCall.put("functions", functions);
        geneCall.put("outsideCall", geneReport.isOutsideCall());
        geneCallList.add(geneCall);
      }
      for (String variant : guideline.getReportVariants()) {
        String call = getGeneReports().stream()
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
        guidelineMap.put("geneCalls", geneCallList);
      }

      guidelineMap.put("reportable", guideline.isReportable());
      guidelineMap.put("uncalledGenes",
          String.join(", ", guideline.getUncalledGenes()));

      guidelineMap.put("matched", guideline.isMatched());
      guidelineMap.put("mutliMatch", guideline.hasMultipleMatches());

      guidelineMap.put("messages", guideline.getMessages().stream()
          .filter(MessageAnnotation.isMessage)
          .map(MessageAnnotation::getMessage)
          .collect(Collectors.toList()));

      guidelineMap.put("footnotes", guideline.getMessages().stream()
          .filter(MessageAnnotation.isFootnote)
          .map(MessageAnnotation::getMessage)
          .collect(Collectors.toList()));

      if (guideline.getCitations().size()>0) {
        guidelineMap.put("citations", guideline.getCitations());
      }

      // special case the display for warfarin recommendation since it's an image
      if (guideline.toString().equals("warfarin")) {
        Map<String,String> imageData = new LinkedHashMap<>();
        imageData.put("url", "http://s3.pgkb.org/attachment/CPIC_warfarin_2017_Fig_2.png");
        imageData.put("altText", "Figure 2 from the CPIC guideline for warfarin");
        guidelineMap.put("image", imageData);
      }

      if (guideline.getMatchingRecommendations() != null) {
        List<Map<String, Object>> groupList = new ArrayList<>();
        for (Recommendation recommendation : guideline.getMatchingRecommendations()) {
          Map<String, Object> groupData = new HashMap<>();

          List<Map<String, String>> annotationList = new ArrayList<>();
          annotationList.add(makeAnnotation("Population", recommendation.getPopulation()));
          for (String geneSymbol : recommendation.getImplications().keySet()) {
            annotationList.add(makeAnnotation("Implication for " + geneSymbol, recommendation.getImplications().get(geneSymbol)));
          }
          for (String geneSymbol : recommendation.getPhenotypes().keySet()) {
            annotationList.add(makeAnnotation("Phenotype for " + geneSymbol, recommendation.getPhenotypes().get(geneSymbol)));
          }
          for (String geneSymbol : recommendation.getActivityScore().keySet()) {
            String score = recommendation.getActivityScore().get(geneSymbol);
            if (score != null && !score.equalsIgnoreCase("n/a")) {
              annotationList.add(makeAnnotation("Activity Score for " + geneSymbol, score));
            }
          }
          for (String geneSymbol : recommendation.getAlleleStatus().keySet()) {
            String status = recommendation.getAlleleStatus().get(geneSymbol);
            if (status != null && !status.equalsIgnoreCase("n/a")) {
              annotationList.add(makeAnnotation("Allele Status for " + geneSymbol, status));
            }
          }
          annotationList.add(makeAnnotation("Recommendation", recommendation.getDrugRecommendation()));
          annotationList.add(makeAnnotation("Classification of Recommendation", recommendation.getClassification()));
          groupData.put("annotations", annotationList);
          groupList.add(groupData);
        }
        guidelineMap.put("groups", groupList);
      }

      sf_logger.debug("Reporting GuidelineReport for {}", guideline.toString());
      guidelines.add(guidelineMap);
    }
    result.put("guidelines", guidelines);


    // Gene calls

    List<Map<String,Object>> geneCallList = new ArrayList<>();
    for (GeneReport geneReport : getGeneReports()) {
      if (geneReport.isIgnored()) {
        continue;
      }

      Map<String,Object> geneCallMap = new HashMap<>();

      geneCallMap.put("gene", geneReport.getGene());

      String phaseStatus;
      if (geneReport.isOutsideCall()) {
        phaseStatus = "Unavailable for calls made outside PharmCAT";
      } else {
        phaseStatus = geneReport.isPhased() ? "Phased" : "Unphased";
      }
      geneCallMap.put("phaseStatus", phaseStatus);

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

  private static Map<String, String> makeAnnotation(String type, String text) {
    Map<String, String> annotationMap = new HashMap<>();
    annotationMap.put("term", type);
    annotationMap.put("annotation", text.replaceAll("[\\n\\r]", " "));
    return annotationMap;
  }
}

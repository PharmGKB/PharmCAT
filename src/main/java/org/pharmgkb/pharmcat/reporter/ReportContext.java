package org.pharmgkb.pharmcat.reporter;

import java.io.IOException;
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
import java.util.function.Predicate;
import java.util.stream.Collectors;
import javax.annotation.Nonnull;
import javax.annotation.Nullable;
import com.google.common.collect.ImmutableList;
import org.apache.commons.lang3.StringUtils;
import org.pharmgkb.pharmcat.definition.IncidentalFinder;
import org.pharmgkb.pharmcat.definition.PhenotypeMap;
import org.pharmgkb.pharmcat.haplotype.DefinitionReader;
import org.pharmgkb.pharmcat.haplotype.model.GeneCall;
import org.pharmgkb.pharmcat.reporter.model.Annotation;
import org.pharmgkb.pharmcat.reporter.model.Group;
import org.pharmgkb.pharmcat.reporter.model.GuidelinePackage;
import org.pharmgkb.pharmcat.reporter.model.MessageAnnotation;
import org.pharmgkb.pharmcat.reporter.model.OutsideCall;
import org.pharmgkb.pharmcat.reporter.model.RelatedGene;
import org.pharmgkb.pharmcat.reporter.model.VariantReport;
import org.pharmgkb.pharmcat.reporter.model.result.GeneReport;
import org.pharmgkb.pharmcat.reporter.model.result.GuidelineReport;
import org.pharmgkb.pharmcat.util.CliUtils;
import org.pharmgkb.pharmcat.util.DataManager;
import org.pharmgkb.pharmcat.util.MessageMatcher;


/**
 * This class acts as a central context for all data needed to generate the final report.
 *
 * It currently gathers
 * <ul>
 *   <li>{@link GeneCall} objects from the named allele matcher</li>
 *   <li>{@link GuidelineReport} objects from dosing guideline annotations</li>
 *   <li>Allele definitions on a per-gene basis</li>
 * </ul>
 *
 * @author greytwist
 * @author Ryan Whaley
 */
public class ReportContext {

  // never display these genes in the gene call list
  private static final List<String> sf_geneBlacklist = ImmutableList.of(
      "G6PD", "HLA-B");
  // never display these drugs in the guideline section
  private static final List<String> sf_drugHidePhenotype = ImmutableList.of(
      "ivacaftor", "peginterferon alfa-2a", "peginterferon alfa-2b", "ribavirin");
  // never display these types of annotations in the guidelines section
  private static final List<String> sf_annotationTermBlacklist = ImmutableList.of(
      "Phenotype (Genotype)", "Metabolizer Status");

  private Map<String,GeneReport> m_geneReports = new TreeMap<>();
  private List<GuidelineReport> m_guidelineReports;
  private PhenotypeMap m_phenotypeMap;
  private IncidentalFinder m_incidentalFinder = new IncidentalFinder();
  private Map<String,String> m_refAlleleForGene = new HashMap<>();

  private final Predicate<String> isGeneIncidental = s -> m_geneReports.values().stream()
      .anyMatch(r -> r.getGene().equals(s) && r.isIncidental());

  /**
   * Public constructor. Compiles all the incoming data into useful objects to be held for later reporting
   * @param calls {@link GeneCall} objects from the sample data
   * @param outsideCalls {@link OutsideCall} objects, non-null but can be empty
   * @param guidelinePackages a List of all the guidelines to try to apply
   */
  public ReportContext(List<GeneCall> calls, @Nonnull List<OutsideCall> outsideCalls, List<GuidelinePackage> guidelinePackages) throws Exception {

    makeGuidelineReports(guidelinePackages);
    makeGeneReports(guidelinePackages);
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
    m_guidelineReports.forEach(r -> r.addMessages(messageMatcher.match(r)));
  }

  /**
   * Takes the raw GuidelinePackage objects from PharmGKG and maps them to {@link GuidelineReport} objects that can be
   * used in the reporter
   * @param guidelinePackages a List of PharmGKB {@link GuidelinePackage} objects
   */
  private void makeGuidelineReports(List<GuidelinePackage> guidelinePackages) {
    m_guidelineReports = guidelinePackages.stream().map(GuidelineReport::new).collect(Collectors.toList());
  }

  /**
   * Makes {@link GeneReport} objects for each of the genes found in PharmGKB {@link GuidelinePackage} objects
   * @param guidelinePackages a List of PharmGKB {@link GuidelinePackage} objects
   */
  private void  makeGeneReports(List<GuidelinePackage> guidelinePackages) {
    for (GuidelinePackage guidelinePackage : guidelinePackages) {
      guidelinePackage.getGuideline().getRelatedGenes().stream()
          .map(RelatedGene::getSymbol)
          .distinct()
          .forEach(s -> m_geneReports.put(s, new GeneReport(s)));
    }
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
          m_incidentalFinder,
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
      GeneReport geneReport = m_geneReports.get(outsideCall.getGene());
      geneReport.setOutsideCallData(outsideCall);

      DiplotypeFactory diplotypeFactory = new DiplotypeFactory(
          outsideCall.getGene(),
          m_phenotypeMap.lookup(outsideCall.getGene()).orElse(null),
          m_incidentalFinder,
          m_refAlleleForGene.get(outsideCall.getGene()));
      geneReport.setDiplotypes(diplotypeFactory, outsideCall);
    }
  }

  /**
   * Adds compiled drug data to each {@link GeneReport} object that's related through a guideline. This gives a useful
   * collection of all drugs that are related to each gene.
   */
  private void compileDrugsForGenes() {
    m_guidelineReports.forEach(r -> {
      for (String gene : r.getRelatedGeneSymbols()) {
        getGeneReport(gene).addRelatedDrugs(r);
      }
    });
  }

  private boolean isReportable(String gene) {
    return m_geneReports.get(gene).isReportable();
  }

  /**
   * Assigns matched guideline groups for all guidelines in this report based on called diplotype functions for each
   * gene and guideline combination.
   */
  private void findMatches() {

    for(GuidelineReport guideline : m_guidelineReports) {
      boolean reportable = guideline.getRelatedGeneSymbols().stream()
          .anyMatch(this::isReportable);

      guideline.setReportable(reportable);

      guideline.getRelatedGeneSymbols().stream()
          .filter(g -> !isReportable(g))
          .forEach(guideline::addUncalledGene);

      if (!reportable) {
        continue;
      }

      makeAllReportGenotypes(guideline);
      guideline.setIncidentalResult(
          guideline.getRelatedGeneSymbols().stream()
              .anyMatch(isGeneIncidental));
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
   * {@link GuidelineReport} and then adds them to the given {@link GuidelineReport}.
   *
   * <em>note:</em> this needs to stay in the {@link ReportContext} since it relies on referencing possibly multiple
   * {@link GeneReport} objects.
   */
  private void makeAllReportGenotypes(GuidelineReport guidelineReport) {
    Set<String> results = new TreeSet<>();
    for (String symbol : guidelineReport.getRelatedGeneSymbols()) {
      results = makeCalledGenotypes(guidelineReport, symbol, results);
    }
    results.forEach(guidelineReport::addReportGenotype);
  }

  /**
   * Makes the genotype Strings for the given gene <code>symbol</code> and either creates a new Set if the passed
   * <code>results</code> is empty or adds to the existing <code>results</code> if there are already entries
   *
   * @param guidelineReport a GuidelineReport with gene calls
   * @param symbol the gene to generate calls for
   * @param results the existing call list to add to
   * @return a new Set of gene calls with the calls for the specified gene
   */
  private Set<String> makeCalledGenotypes(GuidelineReport guidelineReport, String symbol, @Nonnull Set<String> results) {
    if (results.size() == 0) {
      return getGeneReport(symbol).getDiplotypeLookupKeys().stream()
          .map(guidelineReport::translateToPhenotype)
          .collect(Collectors.toSet());
    }
    else {
      Set<String> newResults = new TreeSet<>();
      for (String geno1 : results) {
        getGeneReport(symbol).getDiplotypeLookupKeys().stream().map(guidelineReport::translateToPhenotype).forEach(
            geno2 -> {
              Set<String> genos = new TreeSet<>();
              genos.add(geno1);
              genos.add(geno2);
              newResults.add(String.join(";", genos));
            });
      }
      return newResults;
    }
  }

  public List<GuidelineReport> getGuidelineReports() {
    return m_guidelineReports;
  }

  public Collection<GeneReport> getGeneReports() {
    return m_geneReports.values();
  }

  @Nonnull
  public GeneReport getGeneReport(@Nonnull String geneSymbol) {
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
    for (GeneReport geneReport : getGeneReports()) {
      String symbol = geneReport.getGene();

      // skip any genes on the blacklist
      if (sf_geneBlacklist.contains(symbol)) {
        continue;
      }

      // skip any uncalled genes
      if (!geneReport.isCalled() && geneReport.getReporterDiplotypes().isEmpty()) {
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
    }
    result.put("genotypes", genotypes);
    result.put("totalGenes", genotypes.size());
    result.put("calledGenes", getGeneReports().stream().filter(GeneReport::isCalled).count());

    // Guidelines section
    List<Map<String,Object>> guidelines = new ArrayList<>();
    for (GuidelineReport guideline : new TreeSet<>(getGuidelineReports())) {

      // don't include guidelines that are only on blacklisted genes
      if (guideline.getRelatedGeneSymbols().size() == 1 && sf_geneBlacklist.containsAll(guideline.getRelatedGeneSymbols())) {
        continue;
      }

      Map<String,Object> guidelineMap = new HashMap<>();

      String drugs = String.join(", ", guideline.getRelatedDrugs());
      guidelineMap.put("drugs", drugs);
      guidelineMap.put("summary", guideline.getSummaryHtml());
      guidelineMap.put("url", guideline.getUrl());
      guidelineMap.put("id", guideline.getId());
      String lastModified = guideline.getLastModified() != null
          ? new SimpleDateFormat("YYYY.MM.dd").format(guideline.getLastModified())
          : "not available";
      guidelineMap.put("lastModified", lastModified);

      List<Map<String,Object>> geneCallList = new ArrayList<>();
      for (String gene : guideline.getRelatedGeneSymbols()) {
        GeneReport geneReport = getGeneReport(gene);
        Map<String,Object> geneCall = new LinkedHashMap<>();
        geneCall.put("gene", gene);
        geneCall.put("diplotypes", String.join(", ", geneReport.printDisplayCalls()));
        geneCall.put("showHighlights", !geneReport.getHighlightedVariants().isEmpty());
        geneCall.put("highlightedVariants", geneReport.getHighlightedVariants());
        geneCall.put("outsideCall", geneReport.isOutsideCall());
        geneCallList.add(geneCall);
      }
      for (String variant : guideline.getReportVariants()) {
        String call = getGeneReports().stream()
            .flatMap(g -> g.getVariantReports().stream())
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
      guidelineMap.put("incidental", guideline.isIncidentalResult());

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

      if (guideline.getId().equals("PA166104949")) {
        Map<String,String> imageData = new LinkedHashMap<>();
        imageData.put("url", "http://s3.pgkb.org/attachment/CPIC_warfarin_2017_Fig_2.png");
        imageData.put("altText", "Figure 2 from the CPIC guideline for warfarin");
        guidelineMap.put("image", imageData);
      }

      if (guideline.getMatchingGroups() != null) {
        List<Map<String, Object>> groupList = new ArrayList<>();
        for (Group group : guideline.getMatchingGroups()) {
          Map<String, Object> groupData = new HashMap<>();
          Collection<String> geneFunctions = guideline.getMatchedDiplotypes().get(group.getId());

          List<Map<String, String>> annotationList = new ArrayList<>();
          if (guideline.getRelatedDrugs().stream().noneMatch(sf_drugHidePhenotype::contains)) {
            String alleleFn = geneFunctions.stream()
                .map(f -> f.replace(";", "<br/>"))
                .collect(Collectors.joining("<br/><br/>"));
            annotationList.add(makeAnnotation("Allele Functionality", alleleFn));
            annotationList.add(makeAnnotation("Phenotype", group.getName()));
          }
          for (Annotation ann : group.getAnnotations()) {
            if (sf_annotationTermBlacklist.contains(ann.getType().getTerm())) {
              continue;
            }
            annotationList.add(makeAnnotation(ann.getType().getTerm(), ann.getMarkdown().getHtml()));
          }
          Map<String, String> strengthMap = new HashMap<>();
          strengthMap.put("term", "Classification of Recommendation");
          strengthMap.put("annotation", group.getStrength() != null ? group.getStrength().getTerm() : "N/A");
          annotationList.add(strengthMap);
          groupData.put("annotations", annotationList);
          groupData.put("name", group.getName());
          groupList.add(groupData);
        }
        guidelineMap.put("groups", groupList);
      }

      guidelines.add(guidelineMap);
    }
    result.put("guidelines", guidelines);


    // Gene calls

    List<Map<String,Object>> geneCallList = new ArrayList<>();
    for (GeneReport geneReport : getGeneReports()) {
      if (sf_geneBlacklist.contains(geneReport.getGene())) {
        continue;
      }

      Map<String,Object> geneCallMap = new HashMap<>();

      geneCallMap.put("gene", geneReport.getGene());
      geneCallMap.put("incidental", geneReport.isIncidental());

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

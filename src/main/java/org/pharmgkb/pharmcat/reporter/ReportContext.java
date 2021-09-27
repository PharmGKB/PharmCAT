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
import java.util.Optional;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import org.apache.commons.lang3.StringUtils;
import org.checkerframework.checker.nullness.qual.Nullable;
import org.pharmgkb.pharmcat.definition.MessageList;
import org.pharmgkb.pharmcat.haplotype.model.GeneCall;
import org.pharmgkb.pharmcat.reporter.model.MessageAnnotation;
import org.pharmgkb.pharmcat.reporter.model.VariantReport;
import org.pharmgkb.pharmcat.reporter.model.cpic.Drug;
import org.pharmgkb.pharmcat.reporter.model.cpic.Recommendation;
import org.pharmgkb.pharmcat.reporter.model.result.Diplotype;
import org.pharmgkb.pharmcat.reporter.model.result.DrugReport;
import org.pharmgkb.pharmcat.reporter.model.result.GeneReport;
import org.pharmgkb.pharmcat.util.CliUtils;


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
  private static final Predicate<String> hasValue = (value) -> StringUtils.isNotBlank(value) && !value.equalsIgnoreCase("n/a");

  private final SortedSet<GeneReport> f_geneReports = new TreeSet<>();
  private final List<DrugReport> m_drugReports = new ArrayList<>();

  /**
   * Public constructor. Compiles all the incoming data into useful objects to be held for later reporting
   * @param geneReports {@link GeneReport} objects, non-null but can be empty
   */
  public ReportContext(Collection<GeneReport> geneReports) throws Exception {
    MessageList messageList = new MessageList();

    // add all existing gene data
    f_geneReports.addAll(geneReports);

    // add gene data for genes that have no sample data or outside calls
    DrugCollection drugCollection = new DrugCollection();
    drugCollection.getAllReportableGenes()
        .forEach(this::addGeneReport);

    for (Drug drug : drugCollection.listReportable()) {
      DrugReport drugReport = new DrugReport(drug);
      // set if this drug can be shown in the report
      drugReport.setReportable(drugReport.getRelatedGeneSymbols().stream().anyMatch(this::isReportable));

      // determine which genes don't have calls for this drug
      drugReport.getRelatedGeneSymbols().stream()
          .filter(g -> !isReportable(g))
          .forEach(drugReport::addUncalledGene);

      // add matching recommendations if this drug is reportable
      if (drugReport.isReportable()) {
        List<Map<String,String>> phenoKeys = makePhenotypeKeys(drugReport.getRelatedGeneSymbols());
        for (Map<String,String> phenoKey : phenoKeys) {
          SortedSet<String> diplotypes = findMatchingDiplotypesForPhenotypeKey(phenoKey);
          drugReport.addReportGenotype(phenoKey, diplotypes);
        }
      }
      m_drugReports.add(drugReport);

      // add the inverse relationship to gene reports
      for (String gene : drugReport.getRelatedGeneSymbols()) {
        getGeneReport(gene).addRelatedDrugs(drugReport);
      }
      messageList.match(drugReport, this);

      // add message to drug when a related gene has a *1 allele
      boolean hasStarOne = drugReport.getRelatedGeneSymbols().stream()
          .flatMap((s) -> getGeneReport(s).getReporterDiplotypes().stream())
          .anyMatch((d) -> d.hasAllele("*1"));
      if (hasStarOne) {
        drugReport.addMessage(new MessageAnnotation(
            MessageAnnotation.TYPE_NOTE,
            "The *1 allele assignment is characterized by the absence of variants that are included in the " +
                "underlying allele definitions by either position being reference or missing."
        ));
      }

      // add a message for any gene that has missing data
      drugReport.getRelatedGeneSymbols().stream()
          .filter((s) -> getGeneReport(s).isMissingVariants().equals(GeneReport.YES))
          .map((s) -> "Some position data used to define " + s + " alleles is missing which may change the matched " +
              "genotype. See the gene section for " + s + " for more information.")
          .forEach((m) -> drugReport.addMessage(new MessageAnnotation(MessageAnnotation.TYPE_NOTE, m)));
    }
  }

  private boolean isReportable(String gene) {
    return findGeneReport(gene)
        .map(GeneReport::isReportable)
        .orElse(false);
  }

  public List<DrugReport> getDrugReports() {
    return m_drugReports;
  }

  public Collection<GeneReport> getGeneReports() {
    return f_geneReports;
  }

  /**
   * Find a {@link GeneReport} based on the gene symbol
   * @param geneSymbol a gene symbol
   */
  public Optional<GeneReport> findGeneReport(String geneSymbol) {
    return getGeneReports().stream().filter(r -> r.getGene().equals(geneSymbol)).findFirst();
  }

  /**
   * Finds the {@link GeneReport} for the given gene symbol and will throw a RuntimeException if it's not found,
   * effectively guaranteeing a non-null result
   * @param geneSymbol a gene symbol to find a report for
   * @return a GeneReport object
   * @throws RuntimeException if the desired gene report does not exist
   */
  public GeneReport getGeneReport(String geneSymbol) {
    return findGeneReport(geneSymbol)
        .orElseThrow(() -> new RuntimeException("No gene exists for " + geneSymbol));
  }

  /**
   * Find a {@link DrugReport} for the drug with the given name.
   * @param drugName the name of the drug to find a report for
   * @return an Optional {@link DrugReport}
   */
  public Optional<DrugReport> findDrugReport(String drugName) {
    return m_drugReports.stream()
        .filter(r -> r.getRelatedDrugs().contains(drugName))
        .findFirst();
  }

  /**
   * Gets a {@link DrugReport} for the drug with the given name. Will throw a {@link RuntimeException} if the drug is
   * not found.
   * @param drugName the name of the drug to find
   * @return a non-null {@link DrugReport}
   */
  public DrugReport getDrugReport(String drugName) {
    return findDrugReport(drugName).orElseThrow(() -> new RuntimeException("No drug exists for " + drugName));
  }

  /**
   * Add a "blank" {@link GeneReport} object just based on the gene symbol if a {@link GeneReport} doesn't already exist
   * @param geneSymbol the gene symbol
   */
  private void addGeneReport(String geneSymbol) {
    if (!findGeneReport(geneSymbol).isPresent()) {
      f_geneReports.add(new GeneReport(geneSymbol));
    }
  }

  /**
   * Make a Map of data that will be used in the final report. This map will be serialized and then applied to the
   * handlebars template.
   * @return a Map of data to serialize into JSON
   */
  public Map<String,Object> compile(@Nullable String title) throws IOException {

    Map<String,Object> result = new HashMap<>();
    result.put("generatedOn", new SimpleDateFormat("MMMMM dd, yyyy").format(new Date()));
    result.put("pharmcatVersion", CliUtils.getVersion());

    if (StringUtils.isNotBlank(title)) {
      result.put("title", title);
    }

    // Genotypes section
    List<Map<String,Object>> genotypes = new ArrayList<>();
    int calledGenes = 0;
    for (GeneReport geneReport : getGeneReports()) {
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
      genotype.put("gene", geneReport.getGeneDisplay());
      genotype.put("called", geneReport.isCalled());
      genotype.put("drugs", geneReport.getRelatedDrugs());
      genotype.put("calls", geneReport.printDisplayCalls());
      genotype.put("functions", geneReport.printDisplayFunctions());
      genotype.put("missingVariants", geneReport.isMissingVariants());
      genotype.put("unphased", !geneReport.isOutsideCall() && !geneReport.isPhased());
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

    // Drugs section
    List<Map<String,Object>> drugReports = new ArrayList<>();
    for (DrugReport drugReport : getDrugReports()) {
      if (!result.containsKey("cpicVersion")) {
        result.put("cpicVersion", drugReport.getCpicVersion());
      }
      Map<String,Object> guidelineMap = new HashMap<>();

      String drugs = String.join(", ", drugReport.getRelatedDrugs());
      guidelineMap.put("drugs", drugs);
      guidelineMap.put("url", drugReport.getUrl());
      guidelineMap.put("id", drugReport.getId());
      String lastModified = drugReport.getLastModified() != null
          ? new SimpleDateFormat("yyyy.MM.dd").format(drugReport.getLastModified())
          : "not available";
      guidelineMap.put("lastModified", lastModified);

      List<Map<String,Object>> geneCallList = new ArrayList<>();
      for (String gene : drugReport.getRelatedGeneSymbols()) {
        findGeneReport(gene).ifPresent((geneReport) -> {
          String functions = geneReport.isCalled() ? String.join("; ", geneReport.printDisplayFunctions()) : null;
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

      guidelineMap.put("reportable", drugReport.isReportable());
      guidelineMap.put("uncalledGenes",
          String.join(", ", drugReport.getUncalledGenes()));

      guidelineMap.put("matched", drugReport.isMatched());
      guidelineMap.put("mutliMatch", drugReport.hasMultipleMatches());

      guidelineMap.put("messages", drugReport.getMessages().stream()
          .filter(MessageAnnotation.isMessage)
          .map(MessageAnnotation::getMessage)
          .collect(Collectors.toList()));

      guidelineMap.put("footnotes", drugReport.getMessages().stream()
          .filter(MessageAnnotation.isFootnote)
          .map(MessageAnnotation::getMessage)
          .collect(Collectors.toList()));

      if (drugReport.getCitations().size()>0) {
        guidelineMap.put("citations", drugReport.getCitations());
      }

      // special case the display for warfarin recommendation since it's an image
      if (drugReport.toString().equals("warfarin")) {
        Map<String,String> imageData = new LinkedHashMap<>();
        imageData.put("url", "https://files.cpicpgx.org/images/warfarin/warfarin_recommendation_diagram.png");
        imageData.put("altText", "Figure 2 from the CPIC guideline for warfarin");
        guidelineMap.put("image", imageData);
      }

      if (drugReport.getMatchingRecommendations() != null) {
        List<Map<String, Object>> groupList = new ArrayList<>();
        for (Recommendation recommendation : drugReport.getMatchingRecommendations()) {
          Map<String, Object> groupData = new HashMap<>();

          List<Map<String, String>> annotationList = new ArrayList<>();
          annotationList.add(makeAnnotation("Population", recommendation.getPopulation()));
          for (String geneSymbol : recommendation.getImplications().keySet()) {
            annotationList.add(makeAnnotation("Implication for " + geneSymbol, recommendation.getImplications().get(geneSymbol)));
          }
          annotationList.add(makeAnnotation("Matched Diplotypes", String.join("; ", recommendation.getMatchedDiplotypes())));
          for (String geneSymbol : recommendation.getPhenotypes().keySet()) {
            annotationList.add(makeAnnotation("Phenotype for " + geneSymbol, recommendation.getPhenotypes().get(geneSymbol)));
          }
          for (String geneSymbol : recommendation.getActivityScore().keySet()) {
            String score = recommendation.getActivityScore().get(geneSymbol);
            if (hasValue.test(score)) {
              annotationList.add(makeAnnotation("Activity Score for " + geneSymbol, score));
            }
          }
          for (String geneSymbol : recommendation.getAlleleStatus().keySet()) {
            String status = recommendation.getAlleleStatus().get(geneSymbol);
            if (hasValue.test(status)) {
              annotationList.add(makeAnnotation("Allele Status for " + geneSymbol, status));
            }
          }
          annotationList.add(makeAnnotation("Recommendation", recommendation.getDrugRecommendation()));
          annotationList.add(makeAnnotation("Classification of Recommendation", recommendation.getClassification()));
          if (hasValue.test(recommendation.getComments())) {
            annotationList.add(makeAnnotation("Comments", recommendation.getComments()));
          }
          groupData.put("annotations", annotationList);
          groupList.add(groupData);
        }
        guidelineMap.put("groups", groupList);
      }

      drugReports.add(guidelineMap);
    }
    result.put("guidelines", drugReports);


    // Gene calls

    List<Map<String,Object>> geneCallList = new ArrayList<>();
    for (GeneReport geneReport : getGeneReports()) {
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

  private static Map<String, String> makeAnnotation(String type, String text) {
    Map<String, String> annotationMap = new HashMap<>();
    annotationMap.put("term", type);
    annotationMap.put("annotation", text.replaceAll("[\\n\\r]", " "));
    return annotationMap;
  }

  /**
   * Makes a List of all possible gene combinations for a given set of genes. For example, if GENEA has 2 calls and
   * GENEB has 1 call then the returned List will contain 2 elements.
   * @param geneSymbols genes to make phenotype keys for
   * @return a List of all possible phenotype keys
   */
  private List<Map<String,String>> makePhenotypeKeys(Collection<String> geneSymbols) {
    List<Map<String,String>> keys = new ArrayList<>();
    for (String geneSymbol : geneSymbols) {
      keys = makePhenotypeKeys(geneSymbol, keys);
    }
    return keys;
  }

  private List<Map<String,String>> makePhenotypeKeys(String geneSymbol, List<Map<String,String>> existingList) {
    if (existingList.isEmpty()) {
      streamGeneLookupKeys(geneSymbol).forEach((k) -> existingList.add(makeGeneLookupMap(geneSymbol, k)));
      return existingList;
    }
    else {
      List<Map<String,String>> newList = new ArrayList<>();
      for (Map<String,String> existingKey : existingList) {
        streamGeneLookupKeys(geneSymbol).forEach((k) -> newList.add(copyGeneLookupMap(existingKey, geneSymbol, k)));
      }
      return newList;
    }
  }

  private Stream<String> streamGeneLookupKeys(String gene) {
    Optional<GeneReport> geneReportOpt = findGeneReport(gene);
    if (geneReportOpt.isPresent()) {
      GeneReport geneReport = geneReportOpt.get();
      if (geneReport.getReporterDiplotypes().size() > 0) {
        return geneReport.getReporterDiplotypes().stream().map(Diplotype::getLookupKey);
      }
    }
    return Stream.of("No Result");
  }

  /**
   * Given a map of Gene-&gt;LookupKey values (e.g. CYP2C19 Normal Metabolizer), find the corresponding diplotypes
   * (e.g. CYP2C19:*1/*2) and return them in a set.
   * @param phenotypeKey a Gene-&gt;LookupKey map from the list of diplotypes
   * @return a SortedSet of each gene:diplotype strings corresponding to the lookup keys in the phenotype key
   */
  private SortedSet<String> findMatchingDiplotypesForPhenotypeKey(Map<String,String> phenotypeKey) {
    SortedSet<String> diplotypes = new TreeSet<>();
    for (String gene : phenotypeKey.keySet()) {
      findGeneReport(gene).ifPresent((r) -> r.getReporterDiplotypes().stream()
          .filter((d) -> d.getLookupKey() != null && d.getLookupKey().equals(phenotypeKey.get(gene)))
          .forEach((d) -> diplotypes.add(gene + ":" + d.printDisplay())));
    }
    return diplotypes;
  }

  private Map<String,String> makeGeneLookupMap(String gene, String lookupKey) {
    Map<String,String> newKey = new HashMap<>();
    newKey.put(gene, lookupKey);
    return newKey;
  }

  private Map<String,String> copyGeneLookupMap(Map<String,String> lookupMap, String gene, String lookupKey) {
    Map<String, String> newKey = new HashMap<>(lookupMap);
    newKey.put(gene, lookupKey);
    return newKey;
  }
}

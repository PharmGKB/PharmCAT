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
import javax.annotation.Nonnull;
import org.apache.commons.lang3.StringUtils;
import org.checkerframework.checker.nullness.qual.Nullable;
import org.pharmgkb.pharmcat.haplotype.model.GeneCall;
import org.pharmgkb.pharmcat.reporter.model.MessageAnnotation;
import org.pharmgkb.pharmcat.reporter.model.VariantReport;
import org.pharmgkb.pharmcat.reporter.model.cpic.Drug;
import org.pharmgkb.pharmcat.reporter.model.cpic.Recommendation;
import org.pharmgkb.pharmcat.reporter.model.result.DrugReport;
import org.pharmgkb.pharmcat.reporter.model.result.GeneReport;
import org.pharmgkb.pharmcat.util.CliUtils;
import org.pharmgkb.pharmcat.util.MessageMatcher;


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
    // add all existing gene data
    f_geneReports.addAll(geneReports);

    // add gene data for genes that have no sample data or outside calls
    DrugCollection drugCollection = new DrugCollection();
    drugCollection.getAllReportableGenes()
        .forEach(this::addGeneReport);

    for (Drug drug : drugCollection) {
      DrugReport drugReport = new DrugReport(drug);
      if (!drugReport.isIgnored()) {
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
            drugReport.addReportGenotype(phenoKey);
          }
        }
        m_drugReports.add(drugReport);

        // add the inverse relationship to gene reports
        for (String gene : drugReport.getRelatedGeneSymbols()) {
          getGeneReport(gene).addRelatedDrugs(drugReport);
        }
      }
    }
  }

  /**
   * Applies the given {@link MessageAnnotation} objects to the data in this report.
   * @param messages a List of {@link MessageAnnotation} objects
   */
  public void applyMessage(List<MessageAnnotation> messages) {
    MessageMatcher messageMatcher = new MessageMatcher(messages, this);

    m_drugReports.forEach(r -> r.addMessages(messageMatcher.match(r)));
  }

  private boolean isReportable(String gene) {
    return getGeneReport(gene).isReportable();
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

  @Nonnull
  public GeneReport getGeneReport(String geneSymbol) {
    return findGeneReport(geneSymbol)
        .orElseThrow(() -> new RuntimeException("No gene exists for " + geneSymbol));
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

    // Drugs section
    List<Map<String,Object>> drugReports = new ArrayList<>();
    for (DrugReport drugReport : getDrugReports()) {
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
        imageData.put("url", "https://s3.pgkb.org/attachment/CPIC_warfarin_2017_Fig_2.png");
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

  private List<Map<String,String>> makePhenotypeKeys(Collection<String> geneSymbols) {
    List<Map<String,String>> keys = new ArrayList<>();
    for (String geneSymbol : geneSymbols) {
      keys = makePhenotypeKeys(geneSymbol, keys);
    }
    return keys;
  }

  private List<Map<String,String>> makePhenotypeKeys(String geneSymbol, List<Map<String,String>> existingList) {
    if (existingList.isEmpty()) {
      findGeneReport(geneSymbol).ifPresent((r) ->
          r.getReporterDiplotypes().forEach((d) -> {
            Map<String,String> newKey = new HashMap<>();
            newKey.put(geneSymbol, d.getLookupKey());
            existingList.add(newKey);
          })
      );
      return existingList;
    }
    else {
      List<Map<String,String>> newList = new ArrayList<>();
      findGeneReport(geneSymbol).ifPresent((r) ->
          r.getReporterDiplotypes().forEach((d) -> {
            for (Map<String,String> existingKey : existingList) {
              Map<String, String> newKey = new HashMap<>(existingKey);
              newKey.put(geneSymbol, d.getLookupKey());
              newList.add(newKey);
            }
          })
      );
      return newList;
    }
  }
}

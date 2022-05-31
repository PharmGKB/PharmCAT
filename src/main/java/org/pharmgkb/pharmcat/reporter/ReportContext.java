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
import java.util.stream.Collectors;
import java.util.stream.Stream;
import javax.annotation.Nonnull;
import com.google.common.base.Preconditions;
import org.apache.commons.lang3.StringUtils;
import org.checkerframework.checker.nullness.qual.Nullable;
import org.pharmgkb.pharmcat.definition.MessageList;
import org.pharmgkb.pharmcat.haplotype.model.GeneCall;
import org.pharmgkb.pharmcat.reporter.model.MessageAnnotation;
import org.pharmgkb.pharmcat.reporter.model.VariantReport;
import org.pharmgkb.pharmcat.reporter.model.cpic.Drug;
import org.pharmgkb.pharmcat.reporter.model.result.Diplotype;
import org.pharmgkb.pharmcat.reporter.model.result.DrugReport;
import org.pharmgkb.pharmcat.reporter.model.result.GeneReport;
import org.pharmgkb.pharmcat.reporter.model.result.Genotype;
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
  private final SortedSet<GeneReport> f_geneReports;
  private final SortedSet<DrugReport> m_drugReports = new TreeSet<>();

  /**
   * Public constructor. Compiles all the incoming data into useful objects to be held for later reporting
   * @param geneReports {@link GeneReport} objects, non-null but can be empty
   */
  public ReportContext(Collection<GeneReport> geneReports) throws IOException {
    MessageList messageList = new MessageList();

    // add GeneReports from the Phenotyper
    f_geneReports = new TreeSet<>(geneReports);

    // get CPIC drug data
    DrugCollection drugCollection = new DrugCollection();
    // get DPWG/PharmGKB drug data
    PgkbGuidelineCollection pgkbGuidelineCollection = new PgkbGuidelineCollection();

    // go through all CPIC drugs
    for (Drug drug : drugCollection.listReportable()) {
      DrugReport drugReport = createOrFindDrugReport(drug);
      drugReport.addDrugData(drug, this);

      // add matching recommendations
      List<Genotype> possibleGenotypes = makePossibleGenotypes(drugReport.getRelatedGeneSymbols());
      for (Genotype genotype : possibleGenotypes) {
        drugReport.matchAnnotationsToGenotype(genotype, drug);
      }
      m_drugReports.add(drugReport);

      // add the inverse relationship to gene reports
      for (String gene : drugReport.getRelatedGeneSymbols()) {
        getGeneReport(gene).addRelatedDrugs(drugReport);
      }
    }

    // go through all DPWG-PharmGKB drugs, we iterate this way because one guideline may have multiple chemicals/drugs
    for (String drugName : pgkbGuidelineCollection.getChemicals()) {
      pgkbGuidelineCollection.findGuidelinePackages(drugName).forEach(guidelinePackage -> {
        DrugReport drugReport = createOrFindDrugReport(drugName);

        // add matching groups for possible genotypes
        List<Genotype> possibleGenotypes = makePossibleGenotypes(guidelinePackage.getGenes());
        for (Genotype genotype : possibleGenotypes) {
          guidelinePackage.match(genotype);
        }

        drugReport.addDrugData(guidelinePackage, this);
        m_drugReports.add(drugReport);
        for (String gene : drugReport.getRelatedGeneSymbols()) {
          getGeneReport(gene).addRelatedDrugs(drugReport);
        }
      });
    }

    // now that all reports are generated, apply the applicable messages
    for (DrugReport drugReport : getDrugReports()) {
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
          .filter((s) -> !getGeneReport(s).isOutsideCall() && getGeneReport(s).isMissingVariants().equals(GeneReport.YES))
          .map((s) -> "Some position data used to define " + s + " alleles is missing which may change the matched " +
              "genotype. See the gene section for " + s + " for more information.")
          .forEach((m) -> drugReport.addMessage(new MessageAnnotation(MessageAnnotation.TYPE_NOTE, m)));
    }
  }

  public SortedSet<DrugReport> getDrugReports() {
    return m_drugReports;
  }

  public Collection<GeneReport> getGeneReports() {
    return f_geneReports;
  }

  private DrugReport createOrFindDrugReport(@Nonnull Drug drug) {
    return createOrFindDrugReport(drug.getDrugName());
  }

  private DrugReport createOrFindDrugReport(@Nonnull String drugName) {
    Preconditions.checkNotNull(drugName);
    return getDrugReports().stream()
        .filter(r -> r.getName().equalsIgnoreCase(drugName))
        .findFirst()
        .orElse(new DrugReport(drugName));
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
  private Optional<DrugReport> findDrugReport(String drugName) {
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
    int totalGenes = 0;
    for (GeneReport geneReport : getGeneReports()) {
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
      genotype.put("calls", geneReport.printDisplayCalls());
      genotype.put("functions", geneReport.printDisplayFunctions());
      genotype.put("missingVariants", geneReport.isMissingVariants());
      genotype.put("unphased", !geneReport.isOutsideCall() && !geneReport.isPhased());
      genotype.put("phenotype", geneReport.printDisplayPhenotypes());
      genotype.put("hasMessages", geneReport.getMessages().size()>0);
      genotype.put("outsideCall", geneReport.isOutsideCall());

      genotypes.add(genotype);

      if (geneReport.isReportable()) {
        calledGenes += 1;
      }
    }
    result.put("genotypes", genotypes);
    result.put("totalGenes", totalGenes);
    result.put("calledGenes", calledGenes);

    // Drugs section
    List<Map<String,Object>> drugReports = new ArrayList<>();
    for (DrugReport drugReport : getDrugReports()) {
      Map<String,Object> drugMap = new LinkedHashMap<>();

      drugMap.put("id", drugReport.getId());
      drugMap.put("name", drugReport.getName());
      drugMap.put("urls", drugReport.getUrls());

      List<Map<String,Object>> geneCallList = new ArrayList<>();
      for (String gene : drugReport.getRelatedGeneSymbols()) {
        findGeneReport(gene).ifPresent((geneReport) -> {
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

  private List<Genotype> makePossibleGenotypes(Collection<String> geneSymbols) {
    List<Diplotype> relevantDiplotypes = new ArrayList<>();
    for (String gene : geneSymbols) {
      GeneReport geneReport = getGeneReport(gene);
      relevantDiplotypes.addAll(geneReport.getReporterDiplotypes());
    }
    return Genotype.makeGenotypes(relevantDiplotypes);
  }
}

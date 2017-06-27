package org.pharmgkb.pharmcat.reporter.io;

import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeSet;
import java.util.stream.Collectors;
import com.google.common.collect.ImmutableList;
import org.pharmgkb.pharmcat.reporter.ReportContext;
import org.pharmgkb.pharmcat.reporter.model.Annotation;
import org.pharmgkb.pharmcat.reporter.model.Group;
import org.pharmgkb.pharmcat.reporter.model.MessageAnnotation;
import org.pharmgkb.pharmcat.reporter.model.VariantReport;
import org.pharmgkb.pharmcat.reporter.model.result.GeneReport;
import org.pharmgkb.pharmcat.reporter.model.result.GuidelineReport;


/**
 * This class collects data from a ReportContext and puts it into a format that the final report can use.
 *
 * @author Ryan Whaley
 */
public class ReportData {

  // never display these genes in the gene call list
  private static final List<String> sf_geneBlacklist = ImmutableList.of("G6PD", "HLA-B");
  private static final List<String> sf_drugHidePhenotype = ImmutableList.of(
      "ivacaftor", "peginterferon alfa-2a", "peginterferon alfa-2b", "ribavirin");
  // never display these types of annotations in the guidelines section
  private static final List<String> sf_annotationTermBlacklist = ImmutableList.of("Phenotype (Genotype)", "Metabolizer Status");

  /**
   * Make a Map that can be used in the final handlebars report
   * @param reportContext a populated report context
   * @return a Map of key to data
   */
  public static Map<String,Object> compile(ReportContext reportContext) throws IOException {

    Map<String,Object> result = new HashMap<>();
    result.put("generatedOn", new SimpleDateFormat("MMMMM dd, yyyy").format(new Date()));

    // Genotypes section
    List<Map<String,Object>> genotypes = new ArrayList<>();
    for (GeneReport geneReport : reportContext.getGeneReports()) {
      String symbol = geneReport.getGene();

      // skip any genes on the blacklist
      if (sf_geneBlacklist.contains(symbol)) {
        continue;
      }

      // skip any uncalled genes
      if (!geneReport.isCalled()) {
        continue;
      }

      Map<String,Object> genotype = new HashMap<>();
      genotype.put("gene", symbol);
      genotype.put("called", geneReport.isCalled());
      genotype.put("drugs", geneReport.getRelatedDrugs());
      genotype.put("calls", geneReport.printDisplayCalls());
      genotype.put("functions", geneReport.printDisplayFunctions());
      genotype.put("uncallableAlleles",
            geneReport.getUncalledHaplotypes() != null && geneReport.getUncalledHaplotypes().size() > 0
      );
      genotype.put("phenotype", geneReport.printDisplayPhenotypes());
      genotype.put("hasMessages", geneReport.getMessages().size()>0);
      genotype.put("astrolabe", geneReport.isAstrolabeCall());

      genotypes.add(genotype);
    }
    result.put("genotypes", genotypes);
    result.put("totalGenes", genotypes.size());
    result.put("calledGenes", reportContext.getGeneReports().stream().filter(GeneReport::isCalled).count());

    // Guidelines section
    List<Map<String,Object>> guidelines = new ArrayList<>();
    for (GuidelineReport guideline : new TreeSet<>(reportContext.getGuidelineReports())) {

      // don't include guidelines that are only on blacklisted genes
      if (guideline.getRelatedGeneSymbols().size() == 1 && guideline.getRelatedGeneSymbols().stream().allMatch(sf_geneBlacklist::contains)) {
        continue;
      }

      Map<String,Object> guidelineMap = new HashMap<>();

      String drugs = guideline.getRelatedDrugs().stream().collect(Collectors.joining(", "));
      guidelineMap.put("drugs", drugs);
      guidelineMap.put("summary", guideline.getSummaryHtml());
      guidelineMap.put("url", guideline.getUrl());
      guidelineMap.put("id", guideline.getId());

      List<Map<String,Object>> geneCallList = new ArrayList<>();
      for (String gene : guideline.getRelatedGeneSymbols()) {
        GeneReport geneReport = reportContext.getGeneReport(gene);
        Map<String,Object> geneCall = new LinkedHashMap<>();
        geneCall.put("gene", gene);
        geneCall.put("diplotypes", geneReport.printDisplayCalls().stream()
            .collect(Collectors.joining(", ")));
        geneCall.put("astrolabe", geneReport.isAstrolabeCall());
        geneCallList.add(geneCall);
      }
      if (geneCallList.size() > 0) {
        guidelineMap.put("geneCalls", geneCallList);
      }

      guidelineMap.put("reportable", guideline.isReportable());
      guidelineMap.put("uncalledGenes",
          guideline.getUncalledGenes().stream().collect(Collectors.joining(", ")));

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

          List<Map<String, String>> annotationList = new ArrayList<>();
          if (guideline.getRelatedDrugs().stream().noneMatch(sf_drugHidePhenotype::contains)) {
            annotationList.add(makeAnnotation("Allele Functionality",
                guideline.getMatchedDiplotypes().get(group.getId()).stream()
                    .collect(Collectors.joining(", "))
            ));
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
    for (GeneReport geneReport : reportContext.getGeneReports()) {
      if (sf_geneBlacklist.contains(geneReport.getGene())) {
        continue;
      }

      Map<String,Object> geneCallMap = new HashMap<>();

      geneCallMap.put("gene", geneReport.getGene());
      geneCallMap.put("incidental", geneReport.isIncidental());

      boolean hasUncalledHaplotypes = geneReport.getUncalledHaplotypes() != null && geneReport.getUncalledHaplotypes().size() > 0;
      geneCallMap.put("hasUncalledHaps", hasUncalledHaplotypes);
      if (hasUncalledHaplotypes) {
        geneCallMap.put("uncalledHaps", geneReport.getUncalledHaplotypes().stream().collect(Collectors.joining(", ")));
      }

      geneCallMap.put("diplotypes", geneReport.printDisplayCalls());

      if (geneReport.getMessages() != null && geneReport.getMessages().size() > 0) {
        geneCallMap.put("warnings", geneReport.getMessages().stream().map(MessageAnnotation::getMessage).collect(Collectors.toList()));
      }

      if (geneReport.getVariantReports().size() > 0) {
        geneCallMap.put("variants", geneReport.getVariantReports());
      } else {
        geneCallMap.put("variantsUnspecified", true);
      }

      geneCallMap.put("astrolabe", geneReport.isAstrolabeCall());

      geneCallMap.put("totalMissingVariants",
          geneReport.getVariantReports().stream().filter(VariantReport::isMissing).count());
      geneCallMap.put("totalVariants", geneReport.getVariantReports().size());

      geneCallMap.put("messages", geneReport.getMessages().stream().map(MessageAnnotation::getMessage).collect(Collectors.toList()));

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

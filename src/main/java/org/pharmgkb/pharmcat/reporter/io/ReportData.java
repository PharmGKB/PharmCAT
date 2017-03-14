package org.pharmgkb.pharmcat.reporter.io;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeSet;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import com.google.common.collect.ImmutableList;
import org.pharmgkb.pharmcat.reporter.ReportContext;
import org.pharmgkb.pharmcat.reporter.model.Annotation;
import org.pharmgkb.pharmcat.reporter.model.Group;
import org.pharmgkb.pharmcat.reporter.model.MessageAnnotation;
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
  // never display these types of annotations in the guidelines section
  private static final List<String> sf_annotationTermBlacklist = ImmutableList.of("Phenotype (Genotype)", "Metabolizer Status");

  /**
   * Make a Map that can be used in the final handlebars report
   * @param reportContext a populated report context
   * @return a Map of key to data
   */
  public static Map<String,Object> compile(ReportContext reportContext) throws IOException {

    Map<String,Object> result = new HashMap<>();

    // Genotypes section
    List<Map<String,Object>> genotypes = new ArrayList<>();
    for (GeneReport geneReport : reportContext.getGeneReports()) {
      String symbol = geneReport.getGene();

      if (sf_geneBlacklist.contains(symbol)) {
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
        Map<String,Object> geneCall = new LinkedHashMap<>();
        geneCall.put("gene", gene);
        geneCall.put("diplotypes", Stream.of(gene).flatMap(reportContext.mapGeneToDiplotypes).collect(Collectors.toList()));
        geneCallList.add(geneCall);
      }
      if (geneCallList.size() > 0) {
        guidelineMap.put("geneCalls", geneCallList);
      }

      guidelineMap.put("reportable", guideline.isReportable());
      guidelineMap.put("uncalledGenes",
          guideline.getUncalledGenes().stream().collect(Collectors.joining(", ")));

      guidelineMap.put("matched", guideline.getMatchingGroups() != null);
      guidelineMap.put("mutliMatch", guideline.getMatchingGroups() != null && guideline.getMatchingGroups().size()>1);
      guidelineMap.put("incidental", guideline.isIncidentalResult());

      guidelineMap.put("messages", guideline.getMessages().stream()
          .filter(MessageAnnotation.isMessage)
          .map(MessageAnnotation::getMessage)
          .collect(Collectors.toList()));

      guidelineMap.put("footnotes", guideline.getMessages().stream()
          .filter(MessageAnnotation.isFootnote)
          .map(MessageAnnotation::getMessage)
          .collect(Collectors.toList()));

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
          groupData.put("matchedPhenotypes", guideline.getMatchedDiplotypes().get(group.getId()).stream()
              .collect(Collectors.joining(", ")));

          List<Map<String, String>> annotationList = new ArrayList<>();
          for (Annotation ann : group.getAnnotations()) {
            if (sf_annotationTermBlacklist.contains(ann.getType().getTerm())) {
              continue;
            }

            Map<String, String> annotationMap = new HashMap<>();
            annotationMap.put("term", ann.getType().getTerm());
            annotationMap.put("annotation", ann.getMarkdown().getHtml().replaceAll("[\\n\\r]", " "));
            annotationList.add(annotationMap);
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

      if (geneReport.getVariants().size() > 0) {
        geneCallMap.put("variants", geneReport.getVariants());
      } else {
        geneCallMap.put("variantsUnspecified", true);
      }

      geneCallMap.put("astrolabe", geneReport.isAstrolabeCall());

      int totalVariants = 0;
      if (geneReport.getMatchData() != null && geneReport.getMatchData().getMissingPositions().size()>0) {
        geneCallMap.put("missingVariants", geneReport.getMatchData().getMissingPositions());
        geneCallMap.put("totalMissingVariants", geneReport.getMatchData().getMissingPositions().size());
        totalVariants += geneReport.getMatchData().getMissingPositions().size();
      }
      totalVariants += geneReport.getVariants().size();
      geneCallMap.put("totalVariants", totalVariants);

      geneCallMap.put("messages", geneReport.getMessages().stream().map(MessageAnnotation::getMessage).collect(Collectors.toList()));

      geneCallList.add(geneCallMap);
    }
    result.put("geneCalls", geneCallList);

    return result;
  }
}

package org.pharmgkb.pharmcat.reporter.io;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeSet;
import java.util.stream.Collectors;
import com.google.common.collect.ImmutableList;
import org.pharmgkb.pharmcat.reporter.ReportContext;
import org.pharmgkb.pharmcat.reporter.model.Annotation;
import org.pharmgkb.pharmcat.reporter.model.Group;
import org.pharmgkb.pharmcat.reporter.model.PharmcatException;
import org.pharmgkb.pharmcat.reporter.model.result.GeneReport;
import org.pharmgkb.pharmcat.reporter.model.result.GuidelineReport;


/**
 * This class collects data from a ReportContext and puts it into a format that the final report can use.
 *
 * @author whaleyr
 */
public class ReportData {

  private static final List<String> sf_geneBlacklist = ImmutableList.of("G6PD", "HLA-B");
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
      genotype.put("call", geneReport.getDips().stream().collect(Collectors.joining(", ")));
      genotype.put("functions", geneReport.getFunctions());
      genotype.put("uncallableAlleles",
            geneReport.getUncalledHaplotypes() != null && geneReport.getUncalledHaplotypes().size() > 0
      );
      genotype.put("phenotype", reportContext.makeGenePhenotypes(symbol));

      genotypes.add(genotype);
    }
    result.put("genotypes", genotypes);
    result.put("totalGenes", genotypes.size());
    result.put("calledGenes", reportContext.getGeneReports().stream().filter(GeneReport::isCalled).count());

    // Guidelines section
    List<Map<String,Object>> guidelines = new ArrayList<>();
    for (GuidelineReport guideline : new TreeSet<>(reportContext.getGuidelineResults())) {
      Map<String,Object> guidelineMap = new HashMap<>();

      String drugs = guideline.getRelatedDrugs().stream().collect(Collectors.joining(", "));
      guidelineMap.put("drugs", drugs);
      guidelineMap.put("summary", guideline.getSummaryHtml());
      guidelineMap.put("url", guideline.getUrl());
      guidelineMap.put("id", guideline.getId());

      guidelineMap.put("diplotypes",
          guideline.getRelatedGeneSymbols().stream()
              .flatMap(reportContext.mapGeneToDiplotypes)
              .collect(Collectors.toList()));

      guidelineMap.put("notReportable", !guideline.isReportable());
      guidelineMap.put("uncalledGenes",
          guideline.getUncalledGenes().stream().collect(Collectors.joining(", ")));

      guidelineMap.put("matched", guideline.getMatchingGroups() != null);
      guidelineMap.put("mutliMatch", guideline.getMatchingGroups() != null && guideline.getMatchingGroups().size()>1);

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

      if (geneReport.getUncalledHaplotypes() != null && geneReport.getUncalledHaplotypes().size() > 0) {
        geneCallMap.put("uncalledHaps", geneReport.getUncalledHaplotypes().stream().collect(Collectors.joining(", ")));
      }

      geneCallMap.put("diplotypes", geneReport.getDips());

      if (geneReport.getExceptionList() != null && geneReport.getExceptionList().size() > 0) {
        geneCallMap.put("warnings", geneReport.getExceptionList().stream().map(PharmcatException::getMessage).collect(Collectors.toList()));
      }

      if (geneReport.getVariants().size() > 0) {
        geneCallMap.put("variants", geneReport.getVariants());
      } else {
        geneCallMap.put("variantsUnspecified", true);
      }

      if (geneReport.getMatchData() != null && geneReport.getMatchData().getMissingPositions().size()>0) {
        geneCallMap.put("missingVariants", geneReport.getMatchData().getMissingPositions());
      }
      geneCallList.add(geneCallMap);
    }
    result.put("geneCalls", geneCallList);

    return result;
  }
}

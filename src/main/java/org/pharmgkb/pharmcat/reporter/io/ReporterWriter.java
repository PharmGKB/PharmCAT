package org.pharmgkb.pharmcat.reporter.io;

import java.io.BufferedWriter;
import java.io.IOException;
import java.lang.invoke.MethodHandles;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import org.pharmgkb.pharmcat.reporter.model.Annotation;
import org.pharmgkb.pharmcat.reporter.model.CPICException;
import org.pharmgkb.pharmcat.reporter.model.Group;
import org.pharmgkb.pharmcat.reporter.resultsJSON.GeneReport;
import org.pharmgkb.pharmcat.reporter.resultsJSON.Interaction;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


public class ReporterWriter {
  private static final Logger sf_logger = LoggerFactory.getLogger(MethodHandles.lookup().lookupClass());
  private static final String sf_outputFileName = "reporter.output.md";

  public static void printResults(Path outputDir, List<Interaction> guidelineResultList, Map<String, GeneReport> symbolToGeneReportMap) throws IOException {
    Path reportPath = outputDir.resolve(sf_outputFileName);
    sf_logger.info("Writing report to {}", reportPath);

    try (BufferedWriter writer = Files.newBufferedWriter(reportPath)) {
      writer.write("# PharmCAT Report\n\n\n");

      writer.write("## Haplotype Calls\n\n");
      for (String gene : symbolToGeneReportMap.keySet()) {
        writer.write("### Gene: " + gene + "\n");
        GeneReport geneReport = symbolToGeneReportMap.get(gene);

        if (geneReport.getDips().size() == 1) {
          writer.write("call: `" + geneReport.getDips().iterator().next() + "`\n");
        } else {
          writer.write("__possible calls__:\n\n");
          for (String dip : geneReport.getDips()) {
            writer.write(" * `" + dip + "`\n");
          }
          writer.write("\n");
        }

        writer.write("__warnings__: ");
        if (geneReport.getExceptionList() == null || geneReport.getExceptionList().size() == 0) {
          writer.write("none applicable\n");
        } else {
          writer.write("\n\n");
          for (CPICException exception : geneReport.getExceptionList()) {
            writer.write(" * "+exception.getMessage()+"\n");
          }
        }
        writer.write("\n");
      }

      writer.write("\n\n");

      writer.write("## Guidelines\n\n");

      for (Interaction guideline : guidelineResultList) {
        writer.write("### Guideline: " + guideline.getName() + "\n");
        writer.write("[guideline on PharmGKB]("+guideline.getUrl()+")\n");

        for (Group group : guideline.getMatchingGroups()) {
          writer.write("_");
          writer.write(group.getName());
          writer.write("_ ");
          writer.write(guideline.getMatchedDiplotypes().get(group.getId()).stream()
              .map(d -> "`"+d+"`")
              .collect(Collectors.joining(", ")));
          writer.write("\n\n");

          writer.write("|Type|Annotation|\n");
          writer.write("|---|---|\n");
          for (Annotation ann : group.getAnnotations()) {
            writer.write("|");
            writer.write(ann.getType().getTerm());
            writer.write("|");
            writer.write(ann.getText().replaceAll("[\\n\\r]", " "));
            writer.write("|\n");
          }
          writer.write("\n");
        }
       writer.write("\n\n");
      }
    }
  }
}

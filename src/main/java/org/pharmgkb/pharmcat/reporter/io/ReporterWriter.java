package org.pharmgkb.pharmcat.reporter.io;

import java.io.BufferedWriter;
import java.io.IOException;
import java.lang.invoke.MethodHandles;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import javax.annotation.Nonnull;
import org.pharmgkb.pharmcat.haplotype.model.json.Variant;
import org.pharmgkb.pharmcat.reporter.model.Annotation;
import org.pharmgkb.pharmcat.reporter.model.CPICException;
import org.pharmgkb.pharmcat.reporter.model.Group;
import org.pharmgkb.pharmcat.reporter.resultsJSON.GeneReport;
import org.pharmgkb.pharmcat.reporter.resultsJSON.Interaction;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


/**
 * Write out a report of information in Markdown format
 */
public class ReporterWriter {
  private static final Logger sf_logger = LoggerFactory.getLogger(MethodHandles.lookup().lookupClass());
  private static final String sf_outputFileName = "annotation.report.md";

  private Path m_outputDir;

  public ReporterWriter(@Nonnull Path outputDir) {
    m_outputDir = outputDir;
  }

  public void print(List<Interaction> guidelineResultList, Map<String, GeneReport> symbolToGeneReportMap) throws IOException {
    Path reportPath = m_outputDir.resolve(sf_outputFileName);
    sf_logger.info("Writing report to {}", reportPath);

    try (BufferedWriter writer = Files.newBufferedWriter(reportPath)) {
      writer.write("# PharmCAT Report\n\n\n");

      writer.write("## Haplotype Calls\n\n");
      for (String gene : symbolToGeneReportMap.keySet()) {
        writer.write("### Gene: " + gene + "\n");
        GeneReport geneReport = symbolToGeneReportMap.get(gene);

        if (geneReport.getDips().size() == 1) {
          writer.write("#### Call\n\n");
          writer.write(escapeMd(geneReport.getDips().iterator().next()) + "\n");
        } else {
          writer.write("#### Possible Calls\n\n");
          for (String dip : geneReport.getDips()) {
            writer.write(" * " + escapeMd(dip) + "\n");
          }
          writer.write("\n");
        }

        writer.write("#### Warnings ");
        if (geneReport.getExceptionList() == null || geneReport.getExceptionList().size() == 0) {
          writer.write("none applicable\n");
        } else {
          writer.write("\n\n");
          int i=0;
          for (CPICException exception : geneReport.getExceptionList()) {
            writer.write(++i + ". " + escapeMd(exception.getMessage())+"\n");
          }
        }
        writer.write("\n\n");

        writer.write("#### Calls at Positions\n\n");
        writer.write("| Position | RSID | Call |\n");
        writer.write("| -------- | ---- | ---- |\n");
        for (Variant v : geneReport.getVariants()) {
          String rsidDisplay = v.getRsid() == null ? "None" : v.getRsid();
          writer.write(String.format("|%d|%s|%s|\n", v.getPosition(), rsidDisplay, v.getVcfCall()));
        }
      }

      writer.write("\n\n");

      writer.write("## Guidelines\n\n");

      for (Interaction guideline : guidelineResultList) {
        writer.write("---------------------\n\n");
        writer.write("### " + guideline.getName() + "\n\n");
        writer.write("For more information see the [full guideline on PharmGKB]("+guideline.getUrl()+").\n\n");

        if (guideline.getMatchingGroups() == null) {
          writer.write("_no matching annotations found_\n");
          continue;
        }

        if (guideline.getMatchingGroups().size()>1) {
          writer.write("_Note: More than one call was made for the applicable gene so multiple annotation groups could be shown_\n\n");
        }

        for (Group group : guideline.getMatchingGroups()) {
          writer.write("#### Annotations for ");
          writer.write(guideline.getMatchedDiplotypes().get(group.getId()).stream()
              .map(ReporterWriter::escapeMd)
              .collect(Collectors.joining(", ")));
          writer.write("\n\n");

          writer.write("|Type|Annotation|\n");
          writer.write("|---|---|\n");
          for (Annotation ann : group.getAnnotations()) {
            writer.write("|");
            writer.write(ann.getType().getTerm());
            writer.write("|");
            writer.write(escapeMd(ann.getText().replaceAll("[\\n\\r]", " ")));
            writer.write("|\n");
          }
          if (group.getStrength() != null) {
            writer.write("|Classification of Recommendation|"+group.getStrength().getTerm()+"|\n");
          }
          writer.write("\n");
        }
        writer.write("\n");
      }
    }
  }

  private static String escapeMd(String string) {
    if (string == null) {
      return null;
    }
    return string.replaceAll("\\*", "\\\\*");
  }
}

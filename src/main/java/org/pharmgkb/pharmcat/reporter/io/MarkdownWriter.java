package org.pharmgkb.pharmcat.reporter.io;

import java.io.BufferedWriter;
import java.io.IOException;
import java.lang.invoke.MethodHandles;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Set;
import java.util.TreeSet;
import java.util.stream.Collectors;
import javax.annotation.Nonnull;
import org.pharmgkb.pharmcat.haplotype.model.Variant;
import org.pharmgkb.pharmcat.reporter.DataUnifier;
import org.pharmgkb.pharmcat.reporter.model.Annotation;
import org.pharmgkb.pharmcat.reporter.model.CPICException;
import org.pharmgkb.pharmcat.reporter.model.Group;
import org.pharmgkb.pharmcat.reporter.model.result.GeneReport;
import org.pharmgkb.pharmcat.reporter.model.result.GuidelineReport;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


/**
 * Write out a report of information in Markdown format
 */
public class MarkdownWriter {
  private static final Logger sf_logger = LoggerFactory.getLogger(MethodHandles.lookup().lookupClass());

  private Path m_outputFile;

  /**
   * public constructor
   * @param outputFile path to write report file output to
   */
  public MarkdownWriter(@Nonnull Path outputFile) {
    m_outputFile = outputFile;
  }

  /**
   * Print out a report file based on data found in the {@link DataUnifier}
   * @param dataUnifier a {@link DataUnifier} with all data needed to report
   * @throws IOException
   */
  public void print(DataUnifier dataUnifier) throws IOException {
    sf_logger.info("Writing report to {}", m_outputFile);

    try (BufferedWriter writer = Files.newBufferedWriter(m_outputFile)) {
      writer.write("# PharmCAT Report\n\n\n");

      writer.write("## Haplotype Calls\n\n");
      for (String gene : dataUnifier.getSymbolToGeneReportMap().keySet()) {
        writer.write("### Gene: " + gene + "\n");
        GeneReport geneReport = dataUnifier.getSymbolToGeneReportMap().get(gene);

        writer.write("#### Call\n\n");
        if (geneReport.getDips().size() == 1) {
          writer.write(escapeMd(geneReport.getDips().iterator().next()) + "\n");
        } else {
          writer.write(geneReport.getDips().size()+" possible calls\n\n");
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

      Set<GuidelineReport> guidelines = new TreeSet<>(dataUnifier.getGuidelineResults());
      for (GuidelineReport guideline : guidelines) {
        writer.write("---------------------\n\n");
        writer.write("### " + guideline.getName() + "\n\n");

        writer.write(guideline.getSummaryHtml());
        writer.write("\n\n");

        writer.write("For more information see the [full guideline on PharmGKB]("+guideline.getUrl()+").");
        writer.write("\n\n");

        if (!guideline.isReportable()) {
          writer.write("_gene calls insufficient to filter annotations, missing ");
          String missingGenes = guideline.getRelatedGeneSymbols().stream()
              .filter(s -> !dataUnifier.getSymbolToGeneReportMap().keySet().contains(s))
              .collect(Collectors.joining(", "));
          writer.write(missingGenes);
          writer.write("_\n\n");
          continue;
        }

        if (guideline.getMatchingGroups() == null) {
          writer.write("_related genes present but no matching annotations found, check genotype calls_\n\n");
          continue;
        }

        if (guideline.getMatchingGroups().size()>1) {
          writer.write("_Note: More than one call was made for the applicable gene so multiple annotation groups could be shown_\n\n");
        }

        for (Group group : guideline.getMatchingGroups()) {
          writer.write("#### Annotations for ");
          writer.write(guideline.getMatchedDiplotypes().get(group.getId()).stream()
              .map(MarkdownWriter::escapeMd)
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

  /**
   * Make text safe for markdown output. Currently just escapes asterisks (*) so they can be visible in markdown.
   * @param string a string of text to be output to markdown
   * @return the same string but with certain characters escaped
   */
  private static String escapeMd(String string) {
    if (string == null) {
      return null;
    }
    return string.replaceAll("\\*", "\\\\*");
  }
}

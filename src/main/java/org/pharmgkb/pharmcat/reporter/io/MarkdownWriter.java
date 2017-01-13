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
import org.pharmgkb.pharmcat.definition.model.VariantLocus;
import org.pharmgkb.pharmcat.haplotype.model.Variant;
import org.pharmgkb.pharmcat.reporter.ReportContext;
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
  private static final String sf_variantRowTemplate = "| %d | %s | %s |\n";
  private static final String sf_margin = "\n\n";

  private Path m_outputFile;

  /**
   * public constructor
   * @param outputFile path to write report file output to
   */
  public MarkdownWriter(@Nonnull Path outputFile) {
    m_outputFile = outputFile;
  }

  /**
   * Print out a report file based on data found in the {@link ReportContext}
   * @param reportContext a {@link ReportContext} with all data needed to report
   */
  public void print(ReportContext reportContext) throws IOException {
    sf_logger.info("Writing report to {}", m_outputFile);

    try (BufferedWriter writer = Files.newBufferedWriter(m_outputFile)) {
      writer.write("# PharmCAT Report");
      writer.write(sf_margin);
      writer.write(sf_margin);

      writer.write("## Genotypes");
      writer.write(sf_margin);

      writer.write("| Gene | Drugs | Call |\n");
      writer.write("| ---- | ----- | ---- |\n");
      for (GeneReport geneReport : reportContext.getGeneReports()) {
        writer.write("| ");
        writer.write(geneReport.getGene());
        writer.write(" | ");
        writer.write(geneReport.getRelatedDrugs().stream().collect(Collectors.joining(", ")));
        writer.write(" | ");
        writer.write(geneReport.getDips().stream()
            .map(MarkdownWriter::escapeMd)
            .collect(Collectors.joining(", ")));
        writer.write(" |\n");
      }
      writer.write(sf_margin);
      writer.write(sf_margin);

      writer.write("## Guidelines");
      writer.write(sf_margin);

      Set<GuidelineReport> guidelines = new TreeSet<>(reportContext.getGuidelineResults());
      for (GuidelineReport guideline : guidelines) {
        String drugs = guideline.getRelatedDrugs().stream().collect(Collectors.joining(", "));

        writer.write("### " + drugs);
        writer.write(sf_margin);

        writer.write(guideline.getSummaryHtml());
        writer.write(sf_margin);

        writer.write("For more information see the [full guideline on PharmGKB]("+guideline.getUrl()+").");
        writer.write(sf_margin);

        if (!guideline.isReportable()) {
          writer.write("_gene calls insufficient to filter annotations, missing ");
          String missingGenes = guideline.getUncalledGenes().stream()
              .collect(Collectors.joining(", "));
          writer.write(missingGenes);
          writer.write("_");
          writer.write(sf_margin);
          continue;
        }

        if (guideline.getMatchingGroups() == null) {
          if (guideline.isReportable()) {
            writer.write("_alleles called for all necessary genes but no matching annotations found. check the guideline_");
          }
          else {
            writer.write("_related genes present but no matching annotations found, check genotype calls_");
          }
          writer.write(sf_margin);
          continue;
        }

        if (guideline.getMatchingGroups().size()>1) {
          writer.write("_Note: More than one call was made for the applicable gene so multiple annotation groups could be shown_");
          writer.write(sf_margin);
        }

        for (Group group : guideline.getMatchingGroups()) {
          writer.write("#### Annotations for ");
          writer.write(guideline.getMatchedDiplotypes().get(group.getId()).stream()
              .map(MarkdownWriter::escapeMd)
              .collect(Collectors.joining(", ")));
          writer.write(sf_margin);

          writer.write("|Type|Annotation|\n");
          writer.write("|---|---|\n");
          for (Annotation ann : group.getAnnotations()) {
            writer.write("|");
            writer.write(ann.getType().getTerm());
            writer.write("|");
            writer.write(escapeMd(ann.getMarkdown().getHtml().replaceAll("[\\n\\r]", " ")));
            writer.write("|\n");
          }
          if (group.getStrength() != null) {
            writer.write("|Classification of Recommendation|"+group.getStrength().getTerm()+"|\n");
          }
          writer.write("\n");
        }
        writer.write("\n");
      }

      writer.write("## Gene Call Details");
      writer.write(sf_margin);

      for (GeneReport geneReport : reportContext.getGeneReports()) {
        writer.write("### " + geneReport.getGene());
        writer.write(sf_margin);

        writer.write("#### Matching Allele Call");
        writer.write(sf_margin);

        if (geneReport.getUncalledHaplotypes() != null && geneReport.getUncalledHaplotypes().size() > 0) {
          writer.write("_The following haplotypes were not considered due to missing variant data_: "
              + escapeMd(geneReport.getUncalledHaplotypes().stream().collect(Collectors.joining(", "))));
          writer.write(sf_margin);
        }

        if (geneReport.getDips().size() == 1) {
          writer.write("Diplotype call: " + escapeMd(geneReport.getDips().iterator().next()));
        } else {
          writer.write(geneReport.getDips().size()+" possible diplotype calls");
          writer.write("\n");
          for (String dip : geneReport.getDips()) {
            writer.write(" * " + escapeMd(dip) + "\n");
          }
        }
        writer.write(sf_margin);

        writer.write("#### Warnings ");
        if (geneReport.getExceptionList() == null || geneReport.getExceptionList().size() == 0) {
          writer.write("(none)");
          writer.write(sf_margin);
          writer.write("_no warnings applicable to this set of variant data or allele calls_\n");
        } else {
          writer.write("("+geneReport.getExceptionList().size()+")");
          writer.write(sf_margin);
          int i=0;
          for (CPICException exception : geneReport.getExceptionList()) {
            writer.write(++i + ". " + escapeMd(exception.getMessage())+"\n");
          }
        }
        writer.write(sf_margin);

        writer.write("#### Calls at Positions");
        writer.write(sf_margin);

        if (geneReport.getVariants().size() > 0) {
          writer.write("| Position | RSID | Call |\n");
          writer.write("| -------- | ---- | ---- |\n");
          for (Variant v : geneReport.getVariants()) {
            String rsidDisplay = v.getRsid() == null ? "None" : v.getRsid();
            writer.write(String.format(sf_variantRowTemplate, v.getPosition(), rsidDisplay, v.getVcfCall().replace("|", "\\|")));
          }
        } else {
          writer.write("*no variant data specified*");
        }
        if (geneReport.getMatchData() != null && geneReport.getMatchData().getMissingPositions().size()>0) {
          for (VariantLocus variant : geneReport.getMatchData().getMissingPositions()) {
            writer.write(String.format(sf_variantRowTemplate, variant.getPosition(), variant.getRsid(), "*missing*"));
          }
        }

        writer.write(sf_margin);
      }

      writer.write(sf_margin);

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

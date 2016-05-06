package org.pharmgkb.pharmcat.reporter.io;

import java.io.BufferedWriter;
import java.io.IOException;
import java.lang.invoke.MethodHandles;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import org.pharmgkb.pharmcat.reporter.resultsJSON.Gene;
import org.pharmgkb.pharmcat.reporter.resultsJSON.Interaction;
import org.pharmgkb.pharmcat.reporter.resultsJSON.MultiGeneInteraction;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


public class ReporterWriter {
  private static final Logger sf_logger = LoggerFactory.getLogger(MethodHandles.lookup().lookupClass());
  private static final String sf_outputFileName = "reporter.output.txt";

  public static void printResults(Path outputDir, List<Gene> geneListToReport, List<MultiGeneInteraction> multiInteractionsToReport) throws IOException {
    Path reportPath = outputDir.resolve(sf_outputFileName);
    sf_logger.info("Writing report to {}", reportPath);

    try (BufferedWriter writer = Files.newBufferedWriter(reportPath)) {
      for (Gene geneToWrite : geneListToReport) {
        writer.write("Gene: " + geneToWrite.getGene());
        writer.write("\n");
        writer.write("Called Diplotypes: " + geneToWrite.getGene());
        writer.write("\n");

        int exceptCount = geneToWrite.getExceptionList().size();
        for (int i = 0; i < exceptCount; i++) {
          if (i == 0) {
            writer.write("Exceptions:");
            writer.write("\n");
          }
          writer.write("Rule: " + geneToWrite.getExceptionList().get(i).getRule_name());
          writer.write("\n");
          writer.write("Version: " + geneToWrite.getExceptionList().get(i).getVersion());
          writer.write("\n");
          writer.write("Matches: " + geneToWrite.getExceptionList().get(i).getMatches());
          writer.write("\n");
          writer.write("Exception type: " + geneToWrite.getExceptionList().get(i).getException_type());
          writer.write("\n");
          writer.write("Message: " + geneToWrite.getExceptionList().get(i).getMessage());
          writer.write("\n");
        }

        int interactionCount = geneToWrite.getInteractionList().size();
        for (int i = 0; i < interactionCount; i++) {
          Interaction toWrite = geneToWrite.getInteractionList().get(i);
          writer.write("Name: " + toWrite.getName());
          writer.write("\n");
          writer.write("Source: " + toWrite.getSource());
          writer.write("\n");
          writer.write("SummaryHtml: " + toWrite.getSummaryHtml());
          writer.write("\n");
          writer.write("Html Text: " + toWrite.getTextHtml());
          writer.write("\n");
        }
        writer.write("\n\n");
      }
    }
  }
}

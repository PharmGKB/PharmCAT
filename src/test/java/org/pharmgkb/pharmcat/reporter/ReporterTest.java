package org.pharmgkb.pharmcat.reporter;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import org.junit.jupiter.api.Test;
import org.pharmgkb.common.util.PathUtils;
import org.pharmgkb.pharmcat.TestUtils;
import org.pharmgkb.pharmcat.phenotype.Phenotyper;
import org.pharmgkb.pharmcat.reporter.model.result.DrugReport;
import org.pharmgkb.pharmcat.reporter.model.result.GeneReport;
import org.pharmgkb.pharmcat.reporter.model.result.GuidelineReport;

import static org.junit.jupiter.api.Assertions.*;


/**
 * Test run the Reporter and check generated data
 *
 * @author Ryan Whaley
 */
class ReporterTest {
  private static final String PHENOTYPER_FILE_PATH = "org/pharmgkb/pharmcat/phenotyper_output.json";

  @Test
  void testCypc2c9VariantPassthrough() throws Exception {

    Reporter reporter = new Reporter();
    reporter.analyze(Phenotyper.readGeneReports(PathUtils.getPathToResource(PHENOTYPER_FILE_PATH)));

    // test the CYP2C9 data
    GeneReport geneReport = reporter.getContext().getGeneReport("CYP2C9");
    assertTrue(geneReport.isReportable());
    assertTrue(geneReport.isCalled());
    assertFalse(geneReport.isOutsideCall());
    assertNotNull(geneReport.getVariantReports());
    assertTrue(
        geneReport.getVariantOfInterestReports().stream()
            .anyMatch(r -> r.getDbSnpId() != null && r.getDbSnpId().equals("rs12777823")),
        "Exemption variant not included in gene report"
    );

    // test that messages were applied for a drug
    DrugReport warfarinReport = reporter.getContext().getDrugReports().stream()
        .filter(d -> d.getRelatedDrugs().contains("warfarin")).findFirst()
        .orElseThrow(() -> new RuntimeException("No warfarin drug report found"));
    assertEquals(4, warfarinReport.getMessages().size());

    // test that recommendations were matched
    DrugReport desfluraneReport = reporter.getContext().getDrugReports().stream()
        .filter(d -> d.getRelatedDrugs().contains("desflurane")).findFirst()
        .orElseThrow(() -> new RuntimeException("No desflurane drug report found"));
    assertEquals(1, desfluraneReport.getGuidelines().stream().filter(GuidelineReport::isMatched).count());
  }
  
  @Test
  void testMain() throws Exception {
    Path outputReportPath = TestUtils.createTempFile("ReporterTest", ".html");
    String[] args = new String[]{
        "-p",
        PathUtils.getPathToResource(PHENOTYPER_FILE_PATH).toAbsolutePath().toString(),
        "-f",
        outputReportPath.toString(),
        "-t",
        "example_title"
    };
    Reporter.main(args);

    File outputFile = outputReportPath.toFile();
    assertTrue(outputFile.exists());
    assertFalse(outputFile.isDirectory());

    assertEquals(
        "<!DOCTYPE html>",
        Files.lines(outputReportPath).findFirst().orElseThrow(() -> new RuntimeException("Report is empty"))
    );
    
    outputFile.deleteOnExit();
  }
}

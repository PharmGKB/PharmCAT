package org.pharmgkb.pharmcat;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import org.junit.Test;
import org.pharmgkb.common.util.PathUtils;
import org.pharmgkb.pharmcat.reporter.Reporter;
import org.pharmgkb.pharmcat.reporter.model.result.GeneReport;

import static org.junit.Assert.*;


/**
 * Test run the Reporter and check generated data
 *
 * @author Ryan Whaley
 */
public class ReporterTest {
  
  private static final String CALL_FILE_PATH = "org/pharmgkb/pharmcat/haplotype/cyp2c9/s1s1.vcf";
  private static final String MATCHER_FILE_PATH = "org/pharmgkb/pharmcat/example_matcher.json";
  private static final String OUTPUT_DIR = "ReporterTest";

  @Test
  public void testCypc2c9VariantPassthrough() throws Exception {

    Path vcfFile = PathUtils.getPathToResource(CALL_FILE_PATH);

    Path tempOutDir = Files.createTempDirectory(OUTPUT_DIR);

    PharmCAT pharmcat = new PharmCAT(tempOutDir, null, null);
    pharmcat.execute(vcfFile, null, null);

    Reporter reporter = pharmcat.getReporter();
    GeneReport geneReport = reporter.getContext().getGeneReport("CYP2C9");

    assertNotNull(geneReport);
    assertNotNull(geneReport.getVariantReports());


    assertTrue(
        "Exemption variant not included in gene report",
        geneReport.getVariantOfInterestReports().stream()
            .anyMatch(r -> r.getDbSnpId() != null && r.getDbSnpId().equals("rs12777823"))
    );
  }
  
  @Test
  public void testMain() throws Exception {
    Path outputPath = Files.createTempDirectory(OUTPUT_DIR).resolve("ReporterTest.html");
    String[] args = new String[]{
        "-c",
        PathUtils.getPathToResource(MATCHER_FILE_PATH).toAbsolutePath().toString(),
        "-o",
        outputPath.toString(),
        "-t",
        "example_title"
    };
    Reporter.main(args);
    
    File outputFile = outputPath.toFile();
    assertTrue(outputFile.exists());
    assertFalse(outputFile.isDirectory());
    
    outputFile.deleteOnExit();
  }
}

package org.pharmgkb.pharmcat.reporter;

import java.io.IOException;
import java.net.URISyntaxException;
import java.nio.file.Path;
import org.junit.BeforeClass;
import org.junit.Test;
import org.pharmgkb.common.util.PathUtils;

import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;


/**
 * Test class for running the reporter to get output.
 *
 * @author Ryan Whaley
 */
public class ReporterTest {

  private static Reporter s_reporter;

  @BeforeClass
  public static void before() throws URISyntaxException, IOException {
    Path annotationsDir = PathUtils.getPathToResource("org/pharmgkb/pharmcat/annotations");
    s_reporter = new Reporter(annotationsDir);
  }

  @Test
  public void reporterTest() throws Exception {
    Path callerFile     = PathUtils.getPathToResource("org/pharmgkb/pharmcat/reporter/test.haplotyper.output.json");
    Path outputFile     = callerFile.getParent().resolve("test.haplotyper.output.md");

    s_reporter.analyze(callerFile).printMarkdown(outputFile);
  }

  @Test
  public void bigSampleTest() throws Exception {
    Path callerFile     = PathUtils.getPathToResource("org/pharmgkb/pharmcat/reporter/big.sample.json");
    Path outputFile     = callerFile.getParent().resolve("big.sample.md");

    s_reporter.analyze(callerFile).printMarkdown(outputFile);
  }

  @Test
  public void bigMissingTest() throws Exception {
    Path callerFile     = PathUtils.getPathToResource("org/pharmgkb/pharmcat/reporter/big.sample.missing.2c19.loc.json");
    Path outputFile     = callerFile.getParent().resolve("big.sample.missing.2c19.loc.md");

    Reporter reporter = s_reporter.analyze(callerFile);
    assertNotNull(reporter.getGuidelineReports());

    reporter.printMarkdown(outputFile);

    assertTrue(
        "atazanavir and UGT1A1 guideline not called",
        reporter.getGuidelineReports().stream()
            .anyMatch(r -> r.isReportable() && r.getName().contains("atazanavir") && !r.getMatchingGroups().isEmpty())
    );

    assertTrue(
        "citalopram and CYP2C19 guideline not called",
        reporter.getGuidelineReports().stream()
            .anyMatch(r -> r.isReportable() && r.getName().contains("citalopram") && !r.getMatchingGroups().isEmpty())
    );

    assertTrue(
        "ivacaftor and CFTR guideline called but should not be, we don't have reference annotated yet",
        reporter.getGuidelineReports().stream()
            .anyMatch(r -> r.isReportable() && r.getName().contains("ivacaftor") && r.getMatchingGroups() == null)
    );
  }

  @Test
  public void testCyp2c19Missing() throws Exception {
    Path callerFile     = PathUtils.getPathToResource("org/pharmgkb/pharmcat/reporter/test.cyp2c19.missing.json");
    Path outputFile     = callerFile.getParent().resolve("test.cyp2c19.missing.md");

    s_reporter.analyze(callerFile).printMarkdown(outputFile);
  }
}

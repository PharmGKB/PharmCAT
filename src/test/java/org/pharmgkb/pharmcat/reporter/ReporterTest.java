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
    Path exceptionsPath = PathUtils.getPathToResource("org/pharmgkb/pharmcat/exceptions.tsv");
    s_reporter = new Reporter(annotationsDir, exceptionsPath);
  }

  @Test
  public void combinedTest() throws Exception {
    Path callerFile     = PathUtils.getPathToResource("org/pharmgkb/pharmcat/reporter/combined.json");
    Path outputFile     = callerFile.getParent().resolve("combined.html");

    s_reporter.analyze(callerFile).printHtml(outputFile);

    assertNotNull("Guideline reports don't exist", s_reporter.getGuidelineReports());

    assertTrue(
        "atazanavir and UGT1A1 guideline not called",
        s_reporter.getGuidelineReports().stream()
            .anyMatch(r -> r.isReportable() && r.getName().contains("atazanavir") && !r.getMatchingGroups().isEmpty())
    );

    assertTrue(
        "citalopram and CYP2C19 guideline not called",
        s_reporter.getGuidelineReports().stream()
            .anyMatch(r -> r.isReportable() && r.getName().contains("citalopram") && !r.getMatchingGroups().isEmpty())
    );

    assertTrue(
        "ivacaftor and CFTR guideline called but should not be, we don't have reference annotated yet",
        s_reporter.getGuidelineReports().stream()
            .anyMatch(r -> r.isReportable() && r.getName().contains("ivacaftor") && r.getMatchingGroups() == null)
    );
  }

  @Test
  public void testCyp2c19Missing() throws Exception {
    Path callerFile     = PathUtils.getPathToResource("org/pharmgkb/pharmcat/reporter/test.cyp2c19.missing.json");
    Path outputFile     = callerFile.getParent().resolve("test.cyp2c19.missing.html");

    s_reporter.analyze(callerFile).printHtml(outputFile);
  }
}

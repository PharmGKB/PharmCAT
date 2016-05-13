package org.pharmgkb.pharmcat.reporter;

import java.io.IOException;
import java.net.URISyntaxException;
import java.nio.file.Path;
import org.junit.BeforeClass;
import org.junit.Test;
import org.pharmgkb.pharmcat.TestUtil;


/**
 * Test class for running the reporter to get output.
 *
 * @author Ryan Whaley
 */
public class ReporterTest {

  private static Reporter s_reporter;

  @BeforeClass
  public static void before() throws URISyntaxException, IOException {
    Path annotationsDir = TestUtil.getFile("org/pharmgkb/pharmcat/annotations");
    s_reporter = new Reporter(annotationsDir.toFile());
  }

  @Test
  public void reporterTest() throws Exception {
    Path callerFile     = TestUtil.getFile("org/pharmgkb/pharmcat/reporter/test.haplotyper.output.json");
    Path outputFile     = callerFile.getParent().resolve("test.haplotyper.output.md");

    s_reporter.analyze(callerFile.toFile()).printMarkdown(outputFile);
  }

  @Test
  public void bigSampleTest() throws Exception {
    Path callerFile     = TestUtil.getFile("org/pharmgkb/pharmcat/reporter/big.sample.json");
    Path outputFile     = callerFile.getParent().resolve("big.sample.md");

    s_reporter.analyze(callerFile.toFile()).printMarkdown(outputFile);
  }

  @Test
  public void bigMissingTest() throws Exception {
    Path callerFile     = TestUtil.getFile("org/pharmgkb/pharmcat/reporter/big.sample.missing.2c19.loc.json");
    Path outputFile     = callerFile.getParent().resolve("big.sample.missing.2c19.loc.md");

    s_reporter.analyze(callerFile.toFile()).printMarkdown(outputFile);
  }
}

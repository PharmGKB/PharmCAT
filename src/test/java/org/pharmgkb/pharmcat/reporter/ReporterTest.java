package org.pharmgkb.pharmcat.reporter;

import java.net.URISyntaxException;
import java.nio.file.Path;
import org.junit.Before;
import org.junit.Test;
import org.pharmgkb.pharmcat.TestUtil;


/**
 * Test class for running the reporter to get output.
 *
 * @author Ryan Whaley
 */
public class ReporterTest {

  private Path m_annotationsDir;

  @Before
  public void before() throws URISyntaxException {
    m_annotationsDir = TestUtil.getFile("org/pharmgkb/pharmcat/annotations");
  }

  @Test
  public void reporterTest() throws Exception {
    Path callerFile     = TestUtil.getFile("org/pharmgkb/pharmcat/reporter/test.haplotyper.output.json");
    Path outputFile     = callerFile.getParent().resolve("test.haplotyper.output.md");

    Reporter reporter = new Reporter(
        m_annotationsDir.toFile(),
        callerFile.toFile(),
        outputFile);
    reporter.run();
  }

  @Test
  public void bigSampleTest() throws Exception {
    Path callerFile     = TestUtil.getFile("org/pharmgkb/pharmcat/reporter/big.sample.json");
    Path outputFile     = callerFile.getParent().resolve("big.sample.md");

    Reporter reporter = new Reporter(
        m_annotationsDir.toFile(),
        callerFile.toFile(),
        outputFile);
    reporter.run();
  }

  @Test
  public void bigMissingTest() throws Exception {
    Path callerFile     = TestUtil.getFile("org/pharmgkb/pharmcat/reporter/big.sample.missing.2c19.loc.json");
    Path outputFile     = callerFile.getParent().resolve("big.sample.missing.2c19.loc.md");

    Reporter reporter = new Reporter(
        m_annotationsDir.toFile(),
        callerFile.toFile(),
        outputFile);
    reporter.run();
  }
}

package org.pharmgkb.pharmcat.reporter;

import java.nio.file.Path;
import org.junit.Test;
import org.pharmgkb.pharmcat.TestUtil;


/**
 * @author Ryan Whaley
 */
public class ReporterTest {

  @Test
  public void reporterTest() throws Exception {

    Path annotationsDir = TestUtil.getFile("org/pharmgkb/pharmcat/annotations");
    Path callerFile     = TestUtil.getFile("org/pharmgkb/pharmcat/test.haplotyper.output.json");
    Path outputDir      = annotationsDir.getParent();

    if (!outputDir.toFile().exists()) {
      //noinspection ResultOfMethodCallIgnored
      outputDir.toFile().mkdir();
    }

    Reporter reporter = new Reporter(
        annotationsDir.toFile(),
        callerFile.toFile(),
        outputDir);
    reporter.run();
  }
}

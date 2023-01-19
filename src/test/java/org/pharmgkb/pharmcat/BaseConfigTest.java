package org.pharmgkb.pharmcat;

import java.lang.invoke.MethodHandles;
import java.nio.file.Path;
import org.junit.jupiter.api.Test;
import org.pharmgkb.common.util.CliHelper;
import org.pharmgkb.common.util.PathUtils;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.contains;
import static org.hamcrest.Matchers.containsString;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;


/**
 * This is a JUnit test for {@link BaseConfig}.
 *
 * @author Mark Woon
 */
class BaseConfigTest {

  @Test
  void testSamples() throws Exception {
    CliHelper cliHelper = new CliHelper(MethodHandles.lookup().lookupClass())
        .addOption("s", "sample", "comma-separated list of sample IDs", false, "sample_id")
        .addOption("S", "sample-file", "file containing a list of sample IDs, one sample at a line", false, "file");
    cliHelper.parse(new String[] {
        "-s",
        "S1,S2 , S3 , S4,S5,S1",
    });
    BaseConfig config = new BaseConfig(cliHelper);
    assertThat(config.samples, contains("S1", "S2", "S3", "S4", "S5"));
    assertEquals(5, config.samples.size());

    cliHelper.parse(new String[] {
        "-s",
        "S1,S2 , S3 , S4,S5,S1",
        "-S", "foo.txt"
    });
    assertThrows(ReportableException.class, () -> {
      try {
        new BaseConfig(cliHelper);
      } catch (ReportableException ex) {
        assertEquals(ex.getMessage(), "Cannot specify both -s and -S");
        throw ex;
      }
    });

    Path sampleFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/BaseConfigTest-samples-bad.txt");
    cliHelper.parse(new String[] {
        "-S", sampleFile.toString()
    });
    assertThrows(ReportableException.class, () -> {
      try {
        new BaseConfig(cliHelper);
      } catch (ReportableException ex) {
        assertThat(ex.getMessage(), containsString("Please remove comma"));
        throw ex;
      }
    });

    sampleFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/BaseConfigTest-samples.txt");
    cliHelper.parse(new String[] {
        "-S", sampleFile.toString()
    });
    config = new BaseConfig(cliHelper);
    assertThat(config.samples, contains("S1", "S2", "S3", "S4"));
    assertEquals(4, config.samples.size());
  }
}

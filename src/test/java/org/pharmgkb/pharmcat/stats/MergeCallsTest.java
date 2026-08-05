package org.pharmgkb.pharmcat.stats;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertTrue;
import static uk.org.webcompere.systemstubs.SystemStubs.tapSystemErr;


/**
 * Test cases to make sure {@link MergeCalls} is working as expected.
 *
 * @author Mark Woon
 */
class MergeCallsTest {

  /**
   * Regression test: hitting the {@code -o1d}/{@code -o1f} mutually-exclusive validation must not call
   * {@link System#exit(int)} directly (which would kill the JVM, including the test JVM) - it should go through
   * {@link org.pharmgkb.pharmcat.util.CliUtils#failIfNotTest(String)} like every other CLI validation failure in
   * this codebase.
   */
  @Test
  void mainRejectsMutuallyExclusiveOutputOptionsWithoutExitingJvm() throws Exception {
    String systemErr = tapSystemErr(() -> MergeCalls.main(new String[] {
        "-i", "does-not-matter", "-o1d", "-o1f"
    }));
    assertTrue(systemErr.contains("mutually exclusive"),
        "expected the validation message to be surfaced instead of the process exiting: " + systemErr);
  }
}

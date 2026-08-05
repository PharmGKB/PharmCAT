package org.pharmgkb.pharmcat.stats;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertDoesNotThrow;


/**
 * Test cases to make sure {@link MergeReports} is working as expected.
 *
 * @author Mark Woon
 */
class MergeReportsTest {

  /**
   * Smoke test: missing the required {@code -i} option must go through the same
   * {@link org.pharmgkb.pharmcat.util.CliUtils#failIfNotTest()} path every other CLI entry point uses, instead of
   * just returning. Under test, {@code failIfNotTest()} is a no-op (by design, so tests don't exit the JVM), so this
   * can't observe the exit-code fix itself - only that the code path is reachable and safe under test.
   */
  @Test
  void mainDoesNotThrowWhenRequiredInputDirMissing() {
    assertDoesNotThrow(() -> MergeReports.main(new String[0]));
  }
}

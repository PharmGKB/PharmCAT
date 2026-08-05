package org.pharmgkb.pharmcat.util;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertDoesNotThrow;


/**
 * JUnit test for {@link GeneDrugSummary}.
 *
 * @author Mark Woon
 */
class GeneDrugSummaryTest {

  /**
   * Regression test: missing the required {@code -o} option must go through {@link CliUtils#failIfNotTest()}
   * instead of calling {@link System#exit(int)} directly (which would kill the JVM, including the test JVM).
   */
  @Test
  void mainDoesNotThrowWhenRequiredOutputDirMissing() {
    assertDoesNotThrow(() -> GeneDrugSummary.main(new String[0]));
  }
}

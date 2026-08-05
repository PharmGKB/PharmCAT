package org.pharmgkb.pharmcat.subsetter;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertDoesNotThrow;


/**
 * JUnit test for {@link Subsetter}.
 *
 * @author Mark Woon
 */
class SubsetterTest {

  /**
   * Smoke test: missing the required {@code -i} option must go through {@link
   * org.pharmgkb.pharmcat.util.CliUtils#failIfNotTest()} instead of just returning silently. Under test,
   * {@code failIfNotTest()} (no-arg) is a no-op by design, and {@code Subsetter.main()}'s own catch blocks don't
   * propagate a {@link org.pharmgkb.pharmcat.ReportableException} either, so this can't observe the exit-code fix
   * itself - only that the code path is reachable and safe under test.
   */
  @Test
  void mainDoesNotThrowWhenRequiredInputDirMissing() {
    assertDoesNotThrow(() -> Subsetter.main(new String[0]));
  }

  /**
   * Smoke test: the {@code -pos}/{@code -a} mutual-exclusivity check must go through {@link
   * org.pharmgkb.pharmcat.util.CliUtils#failIfNotTest(String)} instead of just returning silently. Same caveat as
   * above - not observable from within JUnit, since the resulting {@link
   * org.pharmgkb.pharmcat.ReportableException} is caught and printed, not propagated, by {@code Subsetter.main()}.
   */
  @Test
  void mainDoesNotThrowWhenPositionsAndAllelesBothSpecified() {
    assertDoesNotThrow(() -> Subsetter.main(new String[] {
        "-i", "does-not-matter", "-pos", "does-not-matter", "-a", "does-not-matter"
    }));
  }
}

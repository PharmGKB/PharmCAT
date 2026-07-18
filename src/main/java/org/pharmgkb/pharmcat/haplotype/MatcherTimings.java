package org.pharmgkb.pharmcat.haplotype;

import java.util.Locale;


/**
 * Small helper for matcher timing diagnostics.
 */
final class MatcherTimings {

  private MatcherTimings() {
  }

  static long start(boolean enabled) {
    if (!enabled) {
      return 0;
    }
    return System.nanoTime();
  }

  static void print(boolean enabled, String context, String stage, long startNanos) {
    if (!enabled) {
      return;
    }
    double elapsedMs = (System.nanoTime() - startNanos) / 1_000_000.0;
    System.out.printf(Locale.ROOT, "NamedAlleleMatcher timing [%s] %s: %.3f ms%n", context, stage, elapsedMs);
  }
}

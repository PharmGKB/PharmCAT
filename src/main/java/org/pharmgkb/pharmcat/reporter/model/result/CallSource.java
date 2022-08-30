package org.pharmgkb.pharmcat.reporter.model.result;

/**
 * Possible values for describing where a genotype call originated from
 */
public enum CallSource {
  NONE,
  MATCHER,
  OUTSIDE;


  /**
   * This comparator prioritizes {@link #OUTSIDE} over {@link #MATCHER}.
   */
  public static int compare(CallSource a, CallSource b) {
    if (a == b) {
      return 0;
    }
    if (a == OUTSIDE) {
      return -1;
    }
    if (b == OUTSIDE) {
      return 1;
    }
    if (a == MATCHER) {
      return -1;
    }
    if (b == MATCHER) {
      return 1;
    }
    if (a == NONE) {
      return -1;
    }
    if (b == NONE) {
      return 1;
    }
    return 0;
  }
}

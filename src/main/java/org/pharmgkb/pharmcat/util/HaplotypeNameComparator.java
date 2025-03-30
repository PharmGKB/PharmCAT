package org.pharmgkb.pharmcat.util;

import java.util.Comparator;
import java.util.List;
import java.util.regex.Pattern;
import org.pharmgkb.pharmcat.haplotype.CombinationMatcher;


/**
 * Comparator for haplotype names that take combination names into account.
 *
 * @author Mark Woon
 */
public class HaplotypeNameComparator implements Comparator<String> {
  public static final Pattern PARTIAL_PATTERN = Pattern.compile("^g\\.\\d+[?=]?[ACGT>\\d]+$");
  private static final Comparator<String> sf_comparator = new HaplotypeNameComparator();

  /**
   * Gets an instance of this comparator.
   *
   * @return an instance of this comparator
   */
  public static Comparator<String> getComparator() {
    return sf_comparator;
  }


  @Override
  public int compare(String name1, String name2) {
    //noinspection StringEquality
    if (name1 == name2) {
      return 0;
    }
    if (name1 == null) {
      return -1;
    }
    if (name2 == null) {
      return 1;
    }
    if (name1.equals(name2)) {
      return 0;
    }

    boolean c1 = CombinationMatcher.isCombinationName(name1);
    boolean c2 = CombinationMatcher.isCombinationName(name2);
    if (c1 != c2) {
      if (c1) {
        return 1;
      }
      return -1;
    }

    if (c1) {
      // both must be combinations
      List<String> h1 = CombinationMatcher.splitCombinationName(name1);
      List<String> h2 = CombinationMatcher.splitCombinationName(name2);

      int rez = Integer.compare(h1.size(), h2.size());
      if (rez != 0) {
        return rez;
      }

      for (int x = 0; x < h1.size(); x += 1) {
        String n1 = h1.get(x);
        String n2 = h2.get(x);

        boolean isPartial1 = PARTIAL_PATTERN.matcher(n1).matches();
        boolean isPartial2 = PARTIAL_PATTERN.matcher(n2).matches();
        if (isPartial1 != isPartial2) {
          // push partials to the bottom
          if (isPartial1) {
            return 1;
          }
          return -1;
        }

        rez = org.pharmgkb.common.comparator.HaplotypeNameComparator.getComparator().compare(n1, n2);
        if (rez != 0) {
          return rez;
        }
      }
      return 0;
    }

    return org.pharmgkb.common.comparator.HaplotypeNameComparator.getComparator().compare(name1, name2);
  }
}

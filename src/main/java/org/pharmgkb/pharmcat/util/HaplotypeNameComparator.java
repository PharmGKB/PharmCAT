package org.pharmgkb.pharmcat.util;

import java.util.Comparator;
import java.util.List;
import org.pharmgkb.pharmcat.haplotype.model.CombinationMatch;

import static org.pharmgkb.pharmcat.haplotype.model.CombinationMatch.COMBINATION_NAME_SPLITTER;
import static org.pharmgkb.pharmcat.haplotype.model.CombinationMatch.extractCombinationName;


/**
 * Comparator for haplotype names that take combination names into account.
 *
 * @author Mark Woon
 */
public class HaplotypeNameComparator implements Comparator<String> {
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

    boolean c1 = CombinationMatch.isCombinationName(name1);
    boolean c2 = CombinationMatch.isCombinationName(name2);
    if (c1 != c2) {
      if (c1) {
        return 1;
      }
      return -1;
    }

    if (c1) {
      List<String> h1 = COMBINATION_NAME_SPLITTER.splitToList(extractCombinationName(name1));
      List<String> h2 = COMBINATION_NAME_SPLITTER.splitToList(extractCombinationName(name2));

      int rez = Integer.compare(h1.size(), h2.size());
      if (rez != 0) {
        return rez;
      }

      for (int x = 0; x < h1.size(); x += 1) {
        String n1 = h1.get(x);
        String n2 = h2.get(x);
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

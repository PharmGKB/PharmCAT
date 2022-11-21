package org.pharmgkb.pharmcat.util;

import java.util.Comparator;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import org.pharmgkb.pharmcat.haplotype.model.CombinationMatch;

import static org.pharmgkb.pharmcat.haplotype.model.CombinationMatch.COMBINATION_NAME_SPLITTER;
import static org.pharmgkb.pharmcat.haplotype.model.CombinationMatch.extractCombinationName;


/**
 * Comparator for haplotype names that take combination names and HGVS names into account.
 * <p>
 * <b>WARNING!  DANGER!</b> This comparator does NOT comply with the comparator contract.
 * It DOES NOT GUARANTEE TRANSITIVITY (i.e. if A > B and B > C, then A > C) if all haplotype names being compared are
 * not of like kind (i.e. all HGVS, all strict star pattern, all loose star pattern, or all non-star pattern).
 * See {@link org.pharmgkb.common.comparator.HaplotypeNameComparator} for details.
 * <p>
 * Support for HGVS names is just a band-aid on this problem to support DPYD.
 *
 * @author Mark Woon
 */
public class HaplotypeNameComparator implements Comparator<String> {
  private static final Comparator<String> sf_comparator = new HaplotypeNameComparator();
  private static final Pattern sf_hgvsPattern = Pattern.compile("(?:.*:)?[cgnmopr]\\.(\\d+).*");

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
    } else if (name2 == null) {
      return 1;
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

      name1 = h1.get(0);
      name2 = h2.get(0);
    }

    // sort HGVS based names nicely
    Matcher h1 = sf_hgvsPattern.matcher(name1);
    Matcher h2 = sf_hgvsPattern.matcher(name2);
    if (h1.matches() && h2.matches()) {
      int p1 = Integer.parseInt(h1.group(1));
      int p2 = Integer.parseInt(h2.group(1));
      if (p1 == p2) {
        return name1.compareTo(name2);
      }
      return Integer.compare(p1, p2);
    }

    return org.pharmgkb.common.comparator.HaplotypeNameComparator.getComparator().compare(name1, name2);
  }
}

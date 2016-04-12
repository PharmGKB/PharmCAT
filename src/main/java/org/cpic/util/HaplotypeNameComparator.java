package org.cpic.util;

/**
 * @author Mark Woon
 */

import org.apache.commons.lang3.ObjectUtils;
import org.apache.commons.lang3.StringUtils;

import java.util.Comparator;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Compare haplotype names, taking star allele nomenclature into account (i.e. treat numbers as numbers, not strings).
 *
 * @author Ryan Whaley
 */
public class HaplotypeNameComparator implements Comparator<String> {
  private static final Pattern sf_starPattern = Pattern.compile("(.*)(\\d+)(.*)");
  private static final Comparator<String> sf_comparator = new HaplotypeNameComparator();

  public static Comparator<String> getComparator() {
    return sf_comparator;
  }

  public int compare(String name1, String name2) {
    Matcher matcher1 = sf_starPattern.matcher(name1);
    Matcher matcher2 = sf_starPattern.matcher(name2);
    //noinspection ResultOfMethodCallIgnored
    matcher1.find();
    //noinspection ResultOfMethodCallIgnored
    matcher2.find();

    if (matcher1.matches() && matcher2.matches()) {
      String prePortion1 = StringUtils.trimToNull(matcher1.group(1));
      String prePortion2 = StringUtils.trimToNull(matcher2.group(1));
      int rez = ObjectUtils.compare(prePortion1, prePortion2);
      if (rez != 0) {
        return rez;
      }

      String starPortion1 = matcher1.group(2);
      String starPortion2 = matcher2.group(2);
      int star1 = Integer.parseInt(starPortion1);
      int star2 = Integer.parseInt(starPortion2);
      rez = ObjectUtils.compare(star1, star2);
      if (rez != 0) {
        return rez;
      }

      String restPortion1 = StringUtils.trimToNull(matcher1.group(3));
      String restPortion2 = StringUtils.trimToNull(matcher2.group(3));
      rez = ObjectUtils.compare(restPortion1, restPortion2);
      if (rez != 0) {
        return rez;
      }
    }
    return ObjectUtils.compare(name1, name2);
  }
}

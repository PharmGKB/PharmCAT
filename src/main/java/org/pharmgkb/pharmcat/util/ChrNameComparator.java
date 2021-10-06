package org.pharmgkb.pharmcat.util;

import java.util.Comparator;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import org.apache.commons.lang3.StringUtils;


/**
 * Comparator for chromosome names in "chrX" format.
 *
 * @author Mark Woon
 */
public class ChrNameComparator implements Comparator<String> {
  public static final ChrNameComparator INSTANCE = new ChrNameComparator();
  private static final Pattern sf_chrPattern = Pattern.compile("chr([0-9XxYyMm]+)");


  @Override
  public int compare(String o1, String o2) {

    Matcher m1 = sf_chrPattern.matcher(o1);
    if (!m1.matches()) {
      throw new IllegalArgumentException(o1 + " is not a valid chr name");
    }
    Matcher m2 = sf_chrPattern.matcher(o2);
    if (!m2.matches()) {
      throw new IllegalArgumentException(o2 + " is not a valid chr name");
    }

    String c1 = m1.group(1);
    String c2 = m2.group(1);

    boolean isName1Numeric = StringUtils.isNumeric(c1);
    boolean isName2Numeric = StringUtils.isNumeric(c2);

    if (isName1Numeric && isName2Numeric) {
      return Integer.compare(Integer.parseInt(c1), Integer.parseInt(c2));
    }
    if (isName1Numeric) {
      return -1;
    }
    if (isName2Numeric) {
      return 1;
    }

    if (c1.equals(c2)) {
      return 0;
    }
    if (c1.equalsIgnoreCase("M")) {
      return 1;
    }
    if (c2.equalsIgnoreCase("M")) {
      return -1;
    }
    return c1.compareTo(c2);
  }
}

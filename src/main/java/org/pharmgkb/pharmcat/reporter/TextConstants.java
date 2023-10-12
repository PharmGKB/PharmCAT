package org.pharmgkb.pharmcat.reporter;

import java.util.List;
import com.google.common.collect.ImmutableList;
import org.apache.commons.lang3.StringUtils;


/**
 * Text values for use in the Reporter.
 *
 * @author Mark Woon
 */
public class TextConstants {
  // n/a needs to be lowercase for JSON map comparisons to work
  public static final String NA = "n/a";
  public static final String SEE_DRUG = "See drug section";
  public static final String UNKNOWN_FUNCTION = "Unknown";
  public static final String UNKNOWN_GENOTYPE = "Empty genotype";
  public static final String NO_RESULT = "No Result";
  public static final String INDETERMINATE = "Indeterminate";
  public static final String UNSPECIFIED = "Unspecified";
  // greater than or equal to symbol
  public static final String GTE = "\u2265";

  /**
   * Displayed when gene has not been called.
   */
  public static final String UNCALLED = "Not called";
  public static final String NO_DATA = UNCALLED + " - no variant data provided";
  public static final List<String> LIST_UNCALLED = ImmutableList.of(UNCALLED);
  public static final List<String> LIST_UNCALLED_NO_DATA = ImmutableList.of(NO_DATA);

  public static final String GENOTYPE_DELIMITER = "/";
  public static final String REFERENCE = "Reference";
  public static final String HOMOZYGOUS_REFERENCE = REFERENCE + GENOTYPE_DELIMITER + REFERENCE;

  /**
   * Detects an unspecified value.
   *
   * @param value the text to test
   * @return true if the text is either blank or {@link #NA}
   */
  public static boolean isUnspecified(String value) {
    return StringUtils.isBlank(value) || StringUtils.stripToEmpty(value).equalsIgnoreCase(NA);
  }
}

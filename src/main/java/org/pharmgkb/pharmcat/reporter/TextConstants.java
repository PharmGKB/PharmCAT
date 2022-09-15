package org.pharmgkb.pharmcat.reporter;

import java.util.List;
import com.google.common.collect.ImmutableList;


/**
 * @author Mark Woon
 */
public class TextConstants {
  public static final String NA = "N/A";
  public static final String SEE_DRUG = "See drug section";
  public static final String UNKNOWN_FUNCTION = "Unknown";

  /**
   * Displayed when gene has not been called.
   */
  public static final String UNCALLED = "Not called";

  public static final List<String> LIST_UNCALLED = ImmutableList.of(UNCALLED);
  public static final List<String> LIST_UNCALLED_NO_DATA = ImmutableList.of(UNCALLED + " - no variant data provided");


}

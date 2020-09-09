package org.pharmgkb.pharmcat.util;

import java.util.regex.Pattern;
import org.apache.commons.lang3.StringUtils;
import org.checkerframework.checker.nullness.qual.Nullable;


/**
 * Code for dealing with variant specific stuff
 *
 * @author Ryan Whaley
 */
public class VariantUtils {

  private static final Pattern sf_validCallPattern = Pattern.compile("^.+[|/].+$");

  /**
   * Check whether a call string is in the expected format.
   * @param call a String representation of a call (can be null)
   * @return true if the call is like "X/X" or "X|X"
   */
  public static boolean isValidCall(@Nullable String call) {
    if (StringUtils.isBlank(call)) {
      return false;
    }
    return sf_validCallPattern.matcher(call).matches();
  }
}

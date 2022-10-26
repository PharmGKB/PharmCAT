package org.pharmgkb.pharmcat.util;

import java.util.Objects;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import org.apache.commons.lang3.StringUtils;
import org.checkerframework.checker.nullness.qual.Nullable;


/**
 * Code for dealing with variant specific stuff.
 *
 * @author Ryan Whaley
 */
public class VariantUtils {

  private static final Pattern sf_validCallPattern = Pattern.compile("^([A-Za-z]+)(?:[|/]([A-Za-z]+))?$");

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

  /**
   * Checks whether the call is heterozygous, e.g. two different alleles.
   * @param call a SNP call
   * @return true if there are two different alleles in the call
   */
  public static boolean isHetCall(@Nullable String call) {
    if (StringUtils.isBlank(call)) {
      return false;
    }
    Matcher m = sf_validCallPattern.matcher(call);
    if (m.find()) {
      if (m.group(2) == null) {
        return false;
      }
      return !Objects.equals(m.group(1), m.group(2));
    }
    return false;
  }
}

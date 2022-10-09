package org.pharmgkb.pharmcat.util;

import java.util.regex.Matcher;
import java.util.regex.Pattern;
import org.apache.commons.lang3.StringUtils;


public class ActivityUtils {
  private static final Pattern decimalPattern = Pattern.compile("[>≥]?\\d+(\\.\\d+)?$");

  /**
   * Normalizes an activity value/score to use decimal notiation if it is appropriate
   * @param value an activity value or score
   * @return a properly formatted activity value or score
   */
  public static String normalize(String value) {
    if (StringUtils.isBlank(value)) {
      return null;
    }
    String trimmed = StringUtils.trim(value);
    Matcher m = decimalPattern.matcher(trimmed);
    if (m.matches()) {
      if (m.group(1) == null) {
        return trimmed + ".0";
      }
    }
    return trimmed;
  }
}

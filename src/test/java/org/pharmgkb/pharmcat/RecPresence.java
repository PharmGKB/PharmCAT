package org.pharmgkb.pharmcat;

import org.apache.commons.lang3.StringUtils;


/**
 * This enumeration is used by tests for indicating whether a recommendation should appear in the PharmCAT report or
 * not.
 *
 * @author Mark Woon
 */
enum RecPresence {
  NO, YES, YES_NO_MATCH;

  static RecPresence fromString(String val) {
    return switch (StringUtils.stripToEmpty(val.toLowerCase())) {
      case "yes" -> YES;
      case "no" -> NO;
      case "no match" -> YES_NO_MATCH;
      default -> throw new IllegalArgumentException("Unsupported value: " + val);
    };
  }
}

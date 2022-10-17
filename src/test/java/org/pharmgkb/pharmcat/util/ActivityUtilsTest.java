package org.pharmgkb.pharmcat.util;

import org.junit.jupiter.api.Test;
import org.pharmgkb.pharmcat.reporter.TextConstants;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNull;


class ActivityUtilsTest {

  @Test
  void testNormalize() {
    assertNull(ActivityUtils.normalize(" "));
    assertNull(ActivityUtils.normalize(null));
    assertEquals("1.0", ActivityUtils.normalize("1.0"));
    assertEquals("1.0", ActivityUtils.normalize("1"));
    assertEquals(TextConstants.GTE + "1.0", ActivityUtils.normalize(TextConstants.GTE + "1"));
    assertEquals(">3.0", ActivityUtils.normalize(">3"));
    assertEquals(TextConstants.NA, ActivityUtils.normalize(TextConstants.NA));
  }
}

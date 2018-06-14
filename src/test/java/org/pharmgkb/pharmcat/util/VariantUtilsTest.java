package org.pharmgkb.pharmcat.util;

import org.junit.Test;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;


/**
 * Test methods in the {@link VariantUtils} class
 *
 * @author Ryan Whaley
 */
public class VariantUtilsTest {

  @Test
  public void testIsValidCall() {
    assertTrue(VariantUtils.isValidCall("A|A"));
    assertTrue(VariantUtils.isValidCall("C/T"));
    assertTrue(VariantUtils.isValidCall("G|TCCCT"));
    assertFalse(VariantUtils.isValidCall("|"));
    assertFalse(VariantUtils.isValidCall("/"));
    assertFalse(VariantUtils.isValidCall(""));
    assertFalse(VariantUtils.isValidCall(null));
  }
}

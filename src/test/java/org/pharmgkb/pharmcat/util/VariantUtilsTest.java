package org.pharmgkb.pharmcat.util;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;


/**
 * Test methods in the {@link VariantUtils} class
 *
 * @author Ryan Whaley
 */
class VariantUtilsTest {

  @Test
  void testIsValidCall() {
    assertTrue(VariantUtils.isValidCall("A|A"));
    assertTrue(VariantUtils.isValidCall("C/T"));
    assertTrue(VariantUtils.isValidCall("G|TCCCT"));
    assertFalse(VariantUtils.isValidCall("|"));
    assertFalse(VariantUtils.isValidCall("/"));
    assertFalse(VariantUtils.isValidCall(""));
    assertFalse(VariantUtils.isValidCall(null));
  }

  @Test
  void testIsHetCall() {
    assertFalse(VariantUtils.isHetCall("A|A"));
    assertFalse(VariantUtils.isHetCall("A/A"));
    assertTrue(VariantUtils.isHetCall("C/T"));
    assertTrue(VariantUtils.isHetCall("C|T"));
    assertTrue(VariantUtils.isHetCall("G|TCCCT"));
    assertFalse(VariantUtils.isHetCall("|"));
    assertFalse(VariantUtils.isHetCall("/"));
    assertFalse(VariantUtils.isHetCall(""));
    assertFalse(VariantUtils.isHetCall(null));
  }
}

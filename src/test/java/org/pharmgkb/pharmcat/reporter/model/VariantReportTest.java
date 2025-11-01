package org.pharmgkb.pharmcat.reporter.model;

import org.junit.jupiter.api.Test;
import org.pharmgkb.pharmcat.haplotype.model.Variant;

import static org.junit.jupiter.api.Assertions.*;


/**
 * Test the {@link VariantReport} object
 *
 * @author Ryan Whaley
 */
class VariantReportTest {

  @Test
  void testSetCall() {
    Variant variant = new Variant(1, "rs0", "A|A", "foo");
    VariantReport variantReport = new VariantReport("chr1", "CFTR", variant);
    assertEquals("A|A", variantReport.getCall());
    assertFalse(variantReport.isMissing());

    Variant secondVariant = new Variant(2, "rs00", "|", "foo");
    VariantReport secondReport = new VariantReport("chr1", "CFTR", secondVariant);
    assertNull(secondReport.getCall());
    assertTrue(secondReport.isMissing());
  }
}

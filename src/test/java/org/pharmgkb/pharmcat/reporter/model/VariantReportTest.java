package org.pharmgkb.pharmcat.reporter.model;

import org.junit.Test;
import org.pharmgkb.pharmcat.haplotype.model.Variant;

import static org.junit.Assert.*;


/**
 * Test the {@link VariantReport} object
 *
 * @author Ryan Whaley
 */
public class VariantReportTest {
  
  @Test
  public void testCallPattern() {
    assertTrue(VariantReport.sf_validCallPattern.matcher("A|A").matches());
    assertTrue(VariantReport.sf_validCallPattern.matcher("C/T").matches());
    assertTrue(VariantReport.sf_validCallPattern.matcher("G|TCCCT").matches());
    assertFalse(VariantReport.sf_validCallPattern.matcher("|").matches());
    assertFalse(VariantReport.sf_validCallPattern.matcher("/").matches());
  }
  
  @Test
  public void testSetCall() {
    Variant variant = new Variant(1, "rs0", "A|A", 1, "foo");
    VariantReport variantReport = new VariantReport("CFTR", variant);
    assertEquals("A|A", variantReport.getCall());
    assertFalse(variantReport.isMissing());

    Variant secondVariant = new Variant(2, "rs00", "|", 2, "foo");
    VariantReport secondReport = new VariantReport("CFTR", secondVariant);
    assertNull(secondReport.getCall());
    assertTrue(secondReport.isMissing());
  }
}

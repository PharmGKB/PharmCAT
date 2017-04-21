package org.pharmgkb.pharmcat.definition;

import java.io.IOException;
import java.util.Collection;
import org.junit.Test;

import static org.junit.Assert.*;


/**
 * Basic check to make sure the {@link VariantAlleleMap} is working.
 *
 * @author Ryan Whaley
 */
public class VariantAlleleMapTest {

  @Test
  public void test() throws IOException {
    VariantAlleleMap variantAlleleMap = new VariantAlleleMap("TPMT");

    Collection<String> alleles = variantAlleleMap.getAlleles(18149045);
    assertNotNull(alleles);
    assertEquals(1, alleles.size());
    assertTrue(alleles.contains("*13"));
  }

  @Test
  public void testFailures() {
    try {
      VariantAlleleMap variantAlleleMap = new VariantAlleleMap("bad-value");
      fail("No definition expected for bad-value");
    }
    catch (IOException ex) {
      // this is fine, everything's fine
    }
  }
}

package org.pharmgkb.pharmcat.reporter.model.result;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotEquals;


/**
 * Test cases to make sure the Haplotype class is working as expected.
 *
 * @author Ryan Whaley
 */
class HaplotypeTest {

  @Test
  void testCyp2c19() {
    Haplotype haplotype = new Haplotype("CYP2C19", "*4A");

    assertEquals("CYP2C19", haplotype.getGene());
    assertEquals("*4A", haplotype.getName());
    assertEquals("*4A", haplotype.printLookup());
    assertEquals("*4A", haplotype.printDisplay());
    assertEquals("CYP2C19*4A", haplotype.toString());
  }

  @Test
  void testCftr() {
    Haplotype haplotype = new Haplotype("CFTR", "ivacaftor non-responsive CFTR sequence");

    assertEquals("CFTR", haplotype.getGene());
    assertEquals("ivacaftor non-responsive CFTR sequence", haplotype.getName());
    assertEquals("ivacaftor non-responsive CFTR sequence", haplotype.printLookup());
    assertEquals("CFTR ivacaftor non-responsive CFTR sequence", haplotype.toString());
  }

  @Test
  void testDpyd() {
    Haplotype haplotype = new Haplotype("DPYD", "Reference");

    assertEquals("DPYD", haplotype.getGene());
    assertEquals("Reference", haplotype.getName());
    assertEquals("Reference", haplotype.printLookup());
    assertEquals("DPYD Reference", haplotype.toString());
  }

  @Test
  void testCftrNonreference() {
    Haplotype haplotype = new Haplotype("CFTR", "D110H");

    assertEquals("CFTR", haplotype.getGene());
    assertEquals("D110H", haplotype.getName());
    assertEquals("D110H", haplotype.printLookup());
    assertEquals("CFTR D110H", haplotype.toString());
  }

  @Test
  void testEquals() {
    Haplotype h1 = new Haplotype("TPMT", "*1");
    Haplotype h2 = new Haplotype("TPMT", "*1");
    Haplotype h3 = new Haplotype("TPMT", "*2");

    assertEquals(h2, h1);
    assertNotEquals(h3, h1);
  }
}

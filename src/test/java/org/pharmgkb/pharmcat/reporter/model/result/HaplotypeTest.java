package org.pharmgkb.pharmcat.reporter.model.result;

import org.junit.jupiter.api.Test;
import org.pharmgkb.pharmcat.definition.IncidentalFinder;

import static org.junit.jupiter.api.Assertions.*;


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
    Haplotype haplotype = new Haplotype("CFTR", "Reference");

    assertEquals("CFTR", haplotype.getGene());
    assertEquals("Reference", haplotype.getName());
    assertEquals("Other", haplotype.printLookup());
    assertEquals("CFTR Reference", haplotype.toString());
  }

  @Test
  void testDpyd() {
    Haplotype haplotype = new Haplotype("DPYD", "Reference");

    assertEquals("DPYD", haplotype.getGene());
    assertEquals("Reference", haplotype.getName());
    assertEquals("Any normal function variant or no variant detected", haplotype.printLookup());
    assertEquals("DPYD Reference", haplotype.toString());
  }

  @Test
  void testIncidental() throws Exception {
    IncidentalFinder incidentalFinder = new IncidentalFinder();

    Haplotype haplotype = new Haplotype("CFTR", "G542X");
    haplotype.setIncidental(incidentalFinder);

    assertEquals("CFTR", haplotype.getGene());
    assertEquals("G542X", haplotype.getName());
    assertEquals("Other", haplotype.printLookup());
    assertEquals("CFTR G542X", haplotype.toString());
    assertTrue(haplotype.isIncidental());
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

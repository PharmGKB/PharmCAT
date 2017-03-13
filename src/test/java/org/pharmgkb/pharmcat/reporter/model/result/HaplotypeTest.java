package org.pharmgkb.pharmcat.reporter.model.result;

import org.junit.Test;
import org.pharmgkb.pharmcat.definition.IncidentalFinder;

import static org.junit.Assert.*;


/**
 * Test cases to make sure the Haplotype class is working as expected.
 *
 * @author Ryan Whaley
 */
public class HaplotypeTest {

  @Test
  public void testCyp2c19() {
    Haplotype haplotype = new Haplotype("CYP2C19", "*4A");

    assertEquals("CYP2C19", haplotype.getGene());
    assertEquals("*4", haplotype.getName());
    assertEquals("*4", haplotype.printLookup());
    assertEquals("CYP2C19*4", haplotype.toString());
  }

  @Test
  public void testCftr() {
    Haplotype haplotype = new Haplotype("CFTR", "Reference");

    assertEquals("CFTR", haplotype.getGene());
    assertEquals("Reference", haplotype.getName());
    assertEquals("Other", haplotype.printLookup());
    assertEquals("CFTR Reference", haplotype.toString());
  }

  @Test
  public void testIncidental() throws Exception {
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
  public void testEquals() {
    Haplotype h1 = new Haplotype("TPMT", "*1");
    Haplotype h2 = new Haplotype("TPMT", "*1");
    Haplotype h3 = new Haplotype("TPMT", "*2");

    assertTrue(h1.equals(h2));
    assertFalse(h1.equals(h3));
  }
}

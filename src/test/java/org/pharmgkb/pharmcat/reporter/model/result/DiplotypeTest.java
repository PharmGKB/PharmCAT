package org.pharmgkb.pharmcat.reporter.model.result;

import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import org.junit.Test;
import org.pharmgkb.pharmcat.definition.IncidentalFinder;

import static org.junit.Assert.*;


/**
 * Test cases to make sure the Diplotype class is working as expected.
 *
 * @author Ryan Whaley
 */
public class DiplotypeTest {

  @Test
  public void testCftr() throws Exception {
    String gene = "CFTR";

    Haplotype h1 = new Haplotype(gene, "Reference");
    h1.setReference(true);
    Haplotype h2 = new Haplotype(gene, "G542X");

    IncidentalFinder incidentalFinder = new IncidentalFinder();
    h1.setIncidental(incidentalFinder);
    h2.setIncidental(incidentalFinder);

    Diplotype diplotype = new Diplotype(gene, h1, h2);

    assertEquals("CFTR:G542X (heterozygous)", diplotype.toString());
    assertEquals("G542X (heterozygous)", diplotype.printBare());
    assertEquals("G542X (heterozygous)", diplotype.printDisplay());
    assertEquals("CFTR:Other/Other", diplotype.printLookupKey());

    Set<String> alleles = diplotype.streamAllelesByZygosity().collect(Collectors.toSet());
    assertEquals(1, alleles.size());
    assertEquals("G542X (heterozygous)", alleles.iterator().next());

    assertTrue(diplotype.hasAllele("G542X"));
    assertTrue(diplotype.hasAllele("Reference"));
    assertTrue(diplotype.isIncidental());
  }

  @Test
  public void testCyp2d6() throws Exception {
    String gene = "CYP2D6";

    Haplotype h1 = new Haplotype(gene, "*1");
    h1.setReference(true);
    Haplotype h2 = new Haplotype(gene, "*3");

    IncidentalFinder incidentalFinder = new IncidentalFinder();
    h1.setIncidental(incidentalFinder);
    h2.setIncidental(incidentalFinder);

    Diplotype diplotype = new Diplotype(gene, h1, h2);

    assertEquals("CYP2D6:*1/*3", diplotype.toString());
    assertEquals("*1/*3", diplotype.printBare());
    assertEquals("*1/*3", diplotype.printDisplay());
    assertEquals("CYP2D6:*1/*3", diplotype.printLookupKey());

    Set<String> alleles = diplotype.streamAllelesByZygosity().collect(Collectors.toSet());
    assertEquals(1, alleles.size());
    assertEquals("*3 (heterozygous)", alleles.iterator().next());

    assertTrue(diplotype.hasAllele("*1"));
    assertTrue(diplotype.hasAllele("*3"));
    assertFalse(diplotype.hasAllele("foo"));
    assertFalse(diplotype.isIncidental());
  }

  @Test
  public void testJoinPhased() {
    String result = Stream.of("*1/*60", "*1/*80")
        .reduce(Diplotype.phasedReducer)
        .orElseThrow(RuntimeException::new);
    assertEquals("*1/*60+*80", result);

    result = Stream.of("*1/*2", "*2/*3")
        .reduce(Diplotype.phasedReducer)
        .orElseThrow(RuntimeException::new);
    assertEquals("*1+*3/*2", result);

    result = Stream.of("*1/*2", "*1/*3", "*2/*3")
        .reduce(Diplotype.phasedReducer)
        .orElseThrow(RuntimeException::new);
    assertEquals("*1+*3/*2+*3", result);

    result = Stream.of("*6/*80+*28", "*60/*80+*28")
        .reduce(Diplotype.phasedReducer)
        .orElseThrow(RuntimeException::new);
    assertEquals("*6+*60/*28+*80", result);
  }
}

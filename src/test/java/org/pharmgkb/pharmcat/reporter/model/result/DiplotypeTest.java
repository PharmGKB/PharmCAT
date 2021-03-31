package org.pharmgkb.pharmcat.reporter.model.result;

import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableMap;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;


/**
 * Test cases to make sure the Diplotype class is working as expected.
 *
 * @author Ryan Whaley
 */
class DiplotypeTest {

  @Test
  void testCftr() {
    String gene = "CFTR";

    Haplotype h1 = new Haplotype(gene, "ivacaftor non-responsive CFTR sequence");
    h1.setReference(true);
    Haplotype h2 = new Haplotype(gene, "D110H");

    Diplotype diplotype = new Diplotype(gene, h1, h2);

    assertEquals("CFTR:D110H (heterozygous)", diplotype.toString());
    assertEquals("D110H (heterozygous)", diplotype.printBare());
    assertEquals("D110H (heterozygous)", diplotype.printDisplay());
    assertEquals(ImmutableMap.of("D110H", 1, "ivacaftor non-responsive CFTR sequence", 1), diplotype.makeLookupMap());

    Set<String> alleles = diplotype.streamAllelesByZygosity().collect(Collectors.toSet());
    assertEquals(1, alleles.size());
    assertEquals("D110H (heterozygous)", alleles.iterator().next());

    assertTrue(diplotype.hasAllele("D110H"));
    assertTrue(diplotype.hasAllele("ivacaftor non-responsive CFTR sequence"));
  }

  @Test
  void testCyp2d6() {
    String gene = "CYP2D6";

    Haplotype h1 = new Haplotype(gene, "*1");
    h1.setReference(true);
    Haplotype h2 = new Haplotype(gene, "*3");

    Diplotype diplotype = new Diplotype(gene, h1, h2);

    assertEquals("CYP2D6:*1/*3", diplotype.toString());
    assertEquals("*1/*3", diplotype.printBare());
    assertEquals("*1/*3", diplotype.printDisplay());
    assertEquals(ImmutableMap.of("*1", 1, "*3", 1), diplotype.makeLookupMap());

    Set<String> alleles = diplotype.streamAllelesByZygosity().collect(Collectors.toSet());
    assertEquals(1, alleles.size());
    assertEquals("*3 (heterozygous)", alleles.iterator().next());

    assertTrue(diplotype.hasAllele("*1"));
    assertTrue(diplotype.hasAllele("*3"));
    assertFalse(diplotype.hasAllele("foo"));
  }

  @Test
  void testJoinPhased() {
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

  @Test
  void testReducePhasedAlleles() {
    String result = Diplotype.reducePhasedDiplotypes(ImmutableList.of("*1/*60", "*1/*80"));
    assertEquals("*1/*60+*80", result);

    result = Diplotype.reducePhasedDiplotypes(ImmutableList.of("*1/*2", "*2/*3"));
    assertEquals("*1+*3/*2", result);

    result = Diplotype.reducePhasedDiplotypes(ImmutableList.of("*1/*2", "*1/*3", "*2/*3")); // throw an error for this scenario
    assertEquals("*1+*3/*2+*3", result);

    result = Diplotype.reducePhasedDiplotypes(ImmutableList.of("*6/*80+*28", "*60/*80+*28"));
    assertEquals("*6+*60/*80+*28", result);

    result = Diplotype.reducePhasedDiplotypes(ImmutableList.of("*6/*80+*28", "*6/*80+*37"));
    assertEquals("*6/(*80+*28)+(*80+*37)", result);

    result = Diplotype.reducePhasedDiplotypes(ImmutableList.of("*6/*80+*28", "*6/*60"));
    assertEquals("*6/*60+(*80+*28)", result);
  }
}

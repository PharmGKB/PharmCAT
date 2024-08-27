package org.pharmgkb.pharmcat.reporter.model.result;

import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableMap;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import org.pharmgkb.pharmcat.DiplotypeUtils;
import org.pharmgkb.pharmcat.Env;
import org.pharmgkb.pharmcat.phenotype.model.OutsideCall;
import org.pharmgkb.pharmcat.reporter.TextConstants;
import org.pharmgkb.pharmcat.reporter.model.DataSource;
import org.pharmgkb.pharmcat.util.HaplotypeNameComparator;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.contains;
import static org.junit.jupiter.api.Assertions.*;


/**
 * Test cases to make sure the Diplotype class is working as expected.
 *
 * @author Ryan Whaley
 */
public class DiplotypeTest {
  private static Env s_env;

  @BeforeAll
  static void prepare() throws Exception {
    s_env = new Env();
  }


  /**
   * Makes a {@link Stream} of zygosity descriptors for a {@link Diplotype}, e.g. *60 (heterozygous). This is a stream
   * since a single Diplotype can be described in 0, 1, or 2 Strings, depending on the particular allele.
   *
   * @return a Stream of 0 or more zygosity Strings
   */
  private Stream<String> streamAllelesByZygosity(Diplotype diplotype) {
    if (diplotype.getAllele1().equals(diplotype.getAllele2())) {
      if (diplotype.getAllele1().isReference()) {
        return Stream.empty();
      }
      return Stream.of(diplotype.getAllele1().getName() + " (homozygous)");
    }
    else {
      Set<String> alleles = new TreeSet<>(HaplotypeNameComparator.getComparator());
      if (!diplotype.getAllele1().isReference()) {
        alleles.add(diplotype.getAllele1().getName() + " (heterozygous)");
      }
      if (diplotype.getAllele2() != null && !diplotype.getAllele2().isReference()) {
        alleles.add(diplotype.getAllele2().getName() + " (heterozygous)");
      }
      return alleles.stream();
    }
  }



  @Test
  void testCftr() {
    String gene = "CFTR";

    Haplotype h1 = new Haplotype(gene, "ivacaftor non-responsive CFTR sequence");
    h1.setReference(true);
    Haplotype h2 = new Haplotype(gene, "D110H");

    Diplotype diplotype = new Diplotype(gene, h1, h2, s_env, DataSource.CPIC);

    assertEquals("CFTR:D110H (heterozygous)", diplotype.toString());
    assertEquals("D110H (heterozygous)", diplotype.getLabel());
    assertEquals(ImmutableMap.of("D110H", 1, "ivacaftor non-responsive CFTR sequence", 1),
        computeLookupMap(diplotype));

    Set<String> alleles = streamAllelesByZygosity(diplotype).collect(Collectors.toSet());
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

    Diplotype diplotype = new Diplotype(gene, h1, h2, s_env, DataSource.CPIC);

    assertEquals("CYP2D6:*1/*3", diplotype.toString());
    assertEquals("*1/*3", diplotype.getLabel());
    assertEquals(ImmutableMap.of("*1", 1, "*3", 1), computeLookupMap(diplotype));

    Set<String> alleles = streamAllelesByZygosity(diplotype).collect(Collectors.toSet());
    assertEquals(1, alleles.size());
    assertEquals("*3 (heterozygous)", alleles.iterator().next());

    assertTrue(diplotype.hasAllele("*1"));
    assertTrue(diplotype.hasAllele("*3"));
    assertFalse(diplotype.hasAllele("foo"));
  }

  @Test
  void testJoinPhased() {
    String result = Stream.of("*1/*60", "*1/*80")
        .reduce(DiplotypeUtils.PhasedReducer)
        .orElseThrow(RuntimeException::new);
    assertEquals("*1/*60+*80", result);

    result = Stream.of("*1/*2", "*2/*3")
        .reduce(DiplotypeUtils.PhasedReducer)
        .orElseThrow(RuntimeException::new);
    assertEquals("*1+*3/*2", result);

    result = Stream.of("*1/*2", "*1/*3", "*2/*3")
        .reduce(DiplotypeUtils.PhasedReducer)
        .orElseThrow(RuntimeException::new);
    assertEquals("*1+*3/*2+*3", result);

    result = Stream.of("*6/*80+*28", "*60/*80+*28")
        .reduce(DiplotypeUtils.PhasedReducer)
        .orElseThrow(RuntimeException::new);
    assertEquals("*6+*60/*28+*80", result);
  }

  @Test
  void testReducePhasedAlleles() {
    String result = DiplotypeUtils.reducePhasedDiplotypes(ImmutableList.of("*1/*60", "*1/*80"));
    assertEquals("*1/*60+*80", result);

    result = DiplotypeUtils.reducePhasedDiplotypes(ImmutableList.of("*1/*2", "*2/*3"));
    assertEquals("*1+*3/*2", result);

    result = DiplotypeUtils.reducePhasedDiplotypes(ImmutableList.of("*1/*2", "*1/*3", "*2/*3")); // throw an error for this scenario
    assertEquals("*1+*3/*2+*3", result);

    result = DiplotypeUtils.reducePhasedDiplotypes(ImmutableList.of("*6/*80+*28", "*60/*80+*28"));
    assertEquals("*6+*60/*80+*28", result);

    result = DiplotypeUtils.reducePhasedDiplotypes(ImmutableList.of("*6/*80+*28", "*6/*80+*37"));
    assertEquals("*6/(*80+*28)+(*80+*37)", result);

    result = DiplotypeUtils.reducePhasedDiplotypes(ImmutableList.of("*6/*80+*28", "*6/*60"));
    assertEquals("*6/*60+(*80+*28)", result);
  }


  @Test
  void testOutsideCall() {
    OutsideCall outsideCall = new OutsideCall(s_env, "CYP2D6\t\tNM\t4.0", 1);
    Diplotype dip = new Diplotype(outsideCall, s_env, DataSource.CPIC);
    //logDiplotype(dip);
    assertTrue(dip.isUnknownAlleles());
    assertThat(dip.getPhenotypes(), contains("Normal Metabolizer"));
    assertEquals("4.0", dip.getActivityScore());
    assertThat(dip.getLookupKeys(), contains("4.0"));
    assertNotNull(dip.getOutsidePhenotypeMismatch());
    assertNotNull(dip.getOutsideActivityScoreMismatch());

    outsideCall = new OutsideCall(s_env, "CYP2D6\t\tNM", 1);
    dip = new Diplotype(outsideCall, s_env, DataSource.CPIC);
    //logDiplotype(dip);
    assertTrue(dip.isUnknownAlleles());
    assertThat(dip.getPhenotypes(), contains("Normal Metabolizer"));
    assertEquals(TextConstants.NA, dip.getActivityScore());
    assertTrue(dip.getLookupKeys().size() > 1);
    assertNull(dip.getOutsidePhenotypeMismatch());
    assertNull(dip.getOutsideActivityScoreMismatch());

    outsideCall = new OutsideCall(s_env, "CYP2D6\t\t\t4.0", 1);
    dip = new Diplotype(outsideCall, s_env, DataSource.CPIC);
    //logDiplotype(dip);
    assertTrue(dip.isUnknownAlleles());
    assertThat(dip.getPhenotypes(), contains("Ultrarapid Metabolizer"));
    assertEquals("4.0", dip.getActivityScore());
    assertThat(dip.getLookupKeys(), contains("4.0"));
    assertNull(dip.getOutsidePhenotypeMismatch());
    assertNull(dip.getOutsideActivityScoreMismatch());


    outsideCall = new OutsideCall(s_env, "CYP2D6\t*1/*1\tNM\t4.0", 1);
    dip = new Diplotype(outsideCall, s_env, DataSource.CPIC);
    //logDiplotype(dip);
    assertNotNull(dip.getAllele1());
    assertEquals("*1", dip.getAllele1().getName());
    assertNotNull(dip.getAllele2());
    assertEquals("*1", dip.getAllele2().getName());
    assertThat(dip.getPhenotypes(), contains("Normal Metabolizer"));
    assertEquals("4.0", dip.getActivityScore());
    assertThat(dip.getLookupKeys(), contains("4.0"));
    assertNotNull(dip.getOutsidePhenotypeMismatch());
    assertNotNull(dip.getOutsideActivityScoreMismatch());

    outsideCall = new OutsideCall(s_env, "CYP2D6\t*1/*1\tNM", 1);
    dip = new Diplotype(outsideCall, s_env, DataSource.CPIC);
    //logDiplotype(dip);
    assertNotNull(dip.getAllele1());
    assertEquals("*1", dip.getAllele1().getName());
    assertNotNull(dip.getAllele2());
    assertEquals("*1", dip.getAllele2().getName());
    assertThat(dip.getPhenotypes(), contains("Normal Metabolizer"));
    assertEquals(TextConstants.NA, dip.getActivityScore());
    assertTrue(dip.getLookupKeys().size() > 1);
    assertNull(dip.getOutsidePhenotypeMismatch());
    assertNull(dip.getOutsideActivityScoreMismatch());


    outsideCall = new OutsideCall(s_env, "CYP2D6\t*1/*1\t\t4.0", 1);
    //logDiplotype(dip);
    dip = new Diplotype(outsideCall, s_env, DataSource.CPIC);
    assertNotNull(dip.getAllele1());
    assertEquals("*1", dip.getAllele1().getName());
    assertNotNull(dip.getAllele2());
    assertEquals("*1", dip.getAllele2().getName());
    assertThat(dip.getPhenotypes(), contains("Ultrarapid Metabolizer"));
    assertEquals("4.0", dip.getActivityScore());
    assertThat(dip.getLookupKeys(), contains("4.0"));
    assertNull(dip.getOutsidePhenotypeMismatch());
    assertNotNull(dip.getOutsideActivityScoreMismatch());
  }

  private void logDiplotype(Diplotype dip) {
    System.out.println();
    System.out.println(dip);
    System.out.println("score: " + dip.getActivityScore());
    System.out.println("score mismatch: " + dip.getOutsideActivityScoreMismatch());
    System.out.println("pheno: " + dip.getPhenotypes());
    System.out.println("pheno mismatch: " + dip.getOutsidePhenotypeMismatch());
    System.out.println("lookup keys: " + dip.getLookupKeys());
  }


  /**
   * Computes the lookup map for this {@link Diplotype}.
   * Meant to be used internally.  This is only public for tests.
   */
  public static Map<String, Integer> computeLookupMap(Diplotype diplotype) {
    return diplotype.computeLookupMap();
  }
}

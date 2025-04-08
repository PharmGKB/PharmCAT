package org.pharmgkb.pharmcat.haplotype.model;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;
import com.google.common.collect.Lists;
import org.junit.jupiter.api.Test;
import org.pharmgkb.pharmcat.definition.model.NamedAllele;
import org.pharmgkb.pharmcat.definition.model.VariantLocus;
import org.pharmgkb.pharmcat.haplotype.MatchData;

import static org.junit.jupiter.api.Assertions.assertEquals;


/**
 * JUnit test for {@link DiplotypeMatch}.
 *
 * @author Mark Woon
 */
class DiplotypeMatchTest {


  @Test
  void testCompareTo() {

    VariantLocus var1 = new VariantLocus("chr1", 1, "g.1T>A");
    VariantLocus var2 = new VariantLocus("chr1", 2, "g.2T>A");
    VariantLocus var3 = new VariantLocus("chr1", 3, "g.3T>A");
    VariantLocus[] variants = new VariantLocus[] { var1, var2, var3 };

    String[] alleles = new String[] { "T", "T", "T" };
    NamedAllele hap1 = new NamedAllele("*1", "*1", alleles, alleles, true);
    hap1.initialize(variants);

    alleles = new String[] { "A", "A", "A" };
    NamedAllele hap2 = new NamedAllele("*4", "*4", alleles, alleles, false);
    hap2.initialize(variants);

    alleles = new String[] { "T", "A", null };
    NamedAllele hap3 = new NamedAllele("*3", "*3", alleles, alleles, false);
    hap3.initialize(variants);

    HaplotypeMatch hm1 = new HaplotypeMatch(hap1);
    HaplotypeMatch hm2 = new HaplotypeMatch(hap2);
    HaplotypeMatch hm3 = new HaplotypeMatch(hap3);

    MatchData dataset = new MatchData("Sample1", "GENE", new TreeMap<>(), variants, null, null, null);

    DiplotypeMatch dm1 = new DiplotypeMatch(hm1, hm1, dataset);
    DiplotypeMatch dm2 = new DiplotypeMatch(hm1, hm2, dataset);
    DiplotypeMatch dm3 = new DiplotypeMatch(hm2, hm2, dataset);
    DiplotypeMatch dm4 = new DiplotypeMatch(hm3, hm2, dataset);

    SortedSet<DiplotypeMatch> matches = new TreeSet<>(Lists.newArrayList(dm1, dm2));
    assertEquals(dm1, matches.first());

    matches = new TreeSet<>(Lists.newArrayList(dm3, dm2));
    assertEquals(dm2, matches.first());

    matches = new TreeSet<>(Lists.newArrayList(dm3, dm4));
    assertEquals(dm3, matches.first());
  }




  /**
   * A - c.1371C>T/[c.557A>G + c.1627A>G (*5) + c.85T>C (*9A)]
   * B - c.1627A>G (*5)/[c.85T>C (*9A) + c.557A>G + c.1371C>T]
   * C - c.85T>C (*9A)/[c.557A>G + c.1371C>T + c.1627A>G (*5)]
   */
  @Test
  void testSorting() {
    NamedAllele na1 = new NamedAllele("na1", "c.1371C>T", new String[0], new String[0], false);
    NamedAllele na2 = new NamedAllele("na2", "c.557A>G", new String[0], new String[0], false);
    NamedAllele na3 = new NamedAllele("na3", "c.1627A>G (*5)", new String[0], new String[0], false);
    NamedAllele na4 = new NamedAllele("na4", "c.85T>C (*9A)", new String[0], new String[0], false);
    MatchData dataset =  new MatchData("Sample1", "GENE", new TreeMap<>(), new VariantLocus[0], null, null, null);

    // c.1371C>T
    HaplotypeMatch hm1 = new HaplotypeMatch(na1);
    // [c.557A>G + c.1627A>G (*5) + c.85T>C (*9A)]
    CombinationMatch cm1 = new CombinationMatch(new VariantLocus[0], "CCCCC", List.of(na3, na4, na2), null);
    DiplotypeMatch dmA = new DiplotypeMatch(hm1, cm1, dataset);

    // c.1627A>G (*5)
    HaplotypeMatch hm2 = new HaplotypeMatch(na3);
    // [c.85T>C (*9A) + c.557A>G + c.1371C>T]
    CombinationMatch cm2 = new CombinationMatch(new VariantLocus[0], "CCCCC", List.of(na4, na2, na1), null);
    DiplotypeMatch dmB = new DiplotypeMatch(hm2, cm2, dataset);

    // c.85T>C (*9A)
    HaplotypeMatch hm3 = new HaplotypeMatch(na4);
    // [c.557A>G + c.1371C>T + c.1627A>G (*5)]
    CombinationMatch cm3 = new CombinationMatch(new VariantLocus[0], "CCCCC", List.of(na2, na1, na3), null);
    DiplotypeMatch dmC = new DiplotypeMatch(hm3, cm3, dataset);

    System.out.println(dmA);
    System.out.println(dmB);
    System.out.println(dmC);
    System.out.println("-------------------");

    List<DiplotypeMatch> matches = new ArrayList<>();
    matches.add(dmA);
    matches.add(dmB);
    matches.add(dmC);
    Collections.sort(matches);
    System.out.println(matches);

    List<DiplotypeMatch> comparisonMatches = new ArrayList<>();
    comparisonMatches.add(dmA);
    comparisonMatches.add(dmC);
    comparisonMatches.add(dmB);
    Collections.sort(comparisonMatches);
    System.out.println(comparisonMatches);
    assertEquals(matches, comparisonMatches);

    comparisonMatches.clear();
    comparisonMatches.add(dmC);
    comparisonMatches.add(dmA);
    comparisonMatches.add(dmB);
    Collections.sort(comparisonMatches);
    System.out.println(comparisonMatches);
    assertEquals(matches, comparisonMatches);

    comparisonMatches.clear();
    comparisonMatches.add(dmC);
    comparisonMatches.add(dmB);
    comparisonMatches.add(dmA);
    Collections.sort(comparisonMatches);
    System.out.println(comparisonMatches);
    assertEquals(matches, comparisonMatches);

    comparisonMatches.clear();
    comparisonMatches.add(dmB);
    comparisonMatches.add(dmA);
    comparisonMatches.add(dmC);
    Collections.sort(comparisonMatches);
    System.out.println(comparisonMatches);
    assertEquals(matches, comparisonMatches);

    comparisonMatches.clear();
    comparisonMatches.add(dmB);
    comparisonMatches.add(dmC);
    comparisonMatches.add(dmA);
    Collections.sort(comparisonMatches);
    System.out.println(comparisonMatches);
    assertEquals(matches, comparisonMatches);
  }
}

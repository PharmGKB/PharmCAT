package org.pharmgkb.pharmcat.haplotype;

import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.TreeMap;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import org.junit.BeforeClass;
import org.junit.Test;
import org.pharmgkb.pharmcat.TestUtil;
import org.pharmgkb.pharmcat.haplotype.model.DiplotypeMatch;
import org.pharmgkb.pharmcat.haplotype.model.HaplotypeMatch;

import static org.junit.Assert.assertEquals;


/**
 * JUnit test for {@link DiplotypeMatcher}.
 *
 * @author Mark Woon
 */
public class DiplotypeMatcherTest {
  private static List<Haplotype> s_haplotypes;


  @BeforeClass
  public static void beforeClass() {

    // initialize test variants
    Variant s_var1 = new Variant("chr1", "test");
    s_var1.setPosition("1");

    Variant s_var2 = new Variant("chr1", "test");
    s_var2.setPosition("2");

    Variant s_var3 = new Variant("chr1", "test");
    s_var3.setPosition("3");

    List<Variant> s_variants = Lists.newArrayList(s_var1, s_var2, s_var3);

    // initialize test haplotypes
    Haplotype s_hap1 = new Haplotype(
        null, "*1", s_variants,
        null,
        Lists.newArrayList("A", "C", "C")
    );
    s_hap1.calculatePermutations(s_variants);

    Haplotype s_hap2 = new Haplotype(
        null, "*4a", Lists.newArrayList(s_var1),
        null,
        Lists.newArrayList("G")
    );
    s_hap2.calculatePermutations(s_variants);

    Haplotype s_hap3 = new Haplotype(
        null, "*4b", s_variants,
        null,
        Lists.newArrayList("G", "T", "T")
    );
    s_hap3.calculatePermutations(s_variants);

    Haplotype s_hap4 = new Haplotype(
        null, "*17", Lists.newArrayList(s_var2, s_var3),
        null,
        Lists.newArrayList("T", "T")
    );
    s_hap4.calculatePermutations(s_variants);

    s_haplotypes = Lists.newArrayList(s_hap1, s_hap2, s_hap3, s_hap4);

    /*
            | 1 | 2 | 3 |
        *1  | A | C | C |
        *4a | G |   |   |
        *4b | G | T | T |
        *17 |   | T | T |
    */
  }


  @Test
  public void test1() throws Exception {

    List<SampleAllele> alleles = Arrays.asList(
        new SampleAllele("chr1", 1, "A", "G", false),
        new SampleAllele("chr1", 2, "C", "C", false),
        new SampleAllele("chr1", 3, "C", "C", false)
    );

    List<DiplotypeMatch> pairMatches = computeHaplotypes(alleles);
    List<String> expectedMatches = Lists.newArrayList("*1/*4a");
    TestUtil.assertDiplotypePairs(expectedMatches, pairMatches);
  }


  @Test
  public void test2() throws Exception {

    List<SampleAllele> alleles = Arrays.asList(
        new SampleAllele("chr1", 1, "A", "G", false),
        new SampleAllele("chr1", 2, "C", "T", false),
        new SampleAllele("chr1", 3, "C", "T", false)
    );

    List<DiplotypeMatch> pairMatches = computeHaplotypes(alleles);
    List<String> expectedMatches = Lists.newArrayList("*1/*4b", "*1/*17", "*1/*4a", "*4a/*17");
    TestUtil.assertDiplotypePairs(expectedMatches, pairMatches);
  }


  @Test
  public void test3() throws Exception {

    List<SampleAllele> alleles = Arrays.asList(
        new SampleAllele("chr1", 1, "A", "A", false),
        new SampleAllele("chr1", 2, "C", "T", false),
        new SampleAllele("chr1", 3, "C", "T", false)
    );

    List<DiplotypeMatch> pairMatches = computeHaplotypes(alleles);
    List<String> expectedMatches = Lists.newArrayList("*1/*17");
    TestUtil.assertDiplotypePairs(expectedMatches, pairMatches);
  }


  @Test
  public void test4() throws Exception {

    List<SampleAllele> alleles = Arrays.asList(
        new SampleAllele("chr1", 1, "G", "G", false),
        new SampleAllele("chr1", 2, "T", "T", false),
        new SampleAllele("chr1", 3, "C", "T", false)
    );

    List<DiplotypeMatch> pairMatches = computeHaplotypes(alleles);
    List<String> expectedMatches = Lists.newArrayList("*4a/*4b", "*4a/*17", "*4a/*4a");
    TestUtil.assertDiplotypePairs(expectedMatches, pairMatches);
  }


  @Test
  public void test5() throws Exception {

    List<SampleAllele> alleles = Arrays.asList(
        new SampleAllele("chr1", 1, "G", "G", false),
        new SampleAllele("chr1", 2, "C", "T", false),
        new SampleAllele("chr1", 3, "C", "T", false)
    );
    List<DiplotypeMatch> pairMatches = computeHaplotypes(alleles);
    List<String> expectedMatches = Lists.newArrayList("*4a/*4b", "*4a/*17", "*4a/*4a");
    TestUtil.assertDiplotypePairs(expectedMatches, pairMatches);
  }


  private List<DiplotypeMatch> computeHaplotypes(List<SampleAllele> alleles) {

    SortedMap<Integer, SampleAllele> sampleAlleleMap = new TreeMap<>();
    for (SampleAllele sa : alleles) {
      sampleAlleleMap.put(sa.getPosition(), sa);
    }
    Set<String> permutations = CombinationUtil.generatePermutations(alleles);

    return new DiplotypeMatcher(sampleAlleleMap, permutations, s_haplotypes).compute();
  }


  @Test
  public void testComparePermutations() throws Exception {

    Set<String> permutations = Sets.newHashSet(
        "1:T;2:A;3:C;4:C;",
        "1:T;2:A;3:C;4:G;",
        "1:T;2:T;3:C;4:C;",
        "1:T;2:T;3:C;4:G;"
    );

    Variant var1 = new Variant("chr1", "test");
    var1.setPosition("1");

    Variant var2 = new Variant("chr1", "test");
    var2.setPosition("2");

    Variant var3 = new Variant("chr1", "test");
    var3.setPosition("3");

    Variant var4 = new Variant("chr1", "test");
    var4.setPosition("4");

    List<Variant> variants = Lists.newArrayList(var1, var2, var3, var4);

    Haplotype hap1 = new Haplotype(
        null, "*1", variants,
        null,
        Lists.newArrayList("T", "A", "C", "C")
    );
    hap1.calculatePermutations(variants);

    Haplotype hap2 = new Haplotype(
        null, "*2", Lists.newArrayList(var2, var3),
        null,
        Lists.newArrayList("T", "C")
    );
    hap2.calculatePermutations(variants);

    Haplotype hap3 = new Haplotype(
        null, "*3", Lists.newArrayList(var4),
        null,
        Lists.newArrayList("GG")
    );
    hap3.calculatePermutations(variants);

    DiplotypeMatcher diplotypeMatcher = new DiplotypeMatcher(new TreeMap<>(), permutations, Lists.newArrayList(hap1, hap2, hap3));

    SortedSet<HaplotypeMatch> matches = diplotypeMatcher.comparePermutations();
    assertEquals(2, matches.size());
    Iterator<HaplotypeMatch> it = matches.iterator();
    assertEquals(hap1, it.next().getHaplotype());
    assertEquals(hap2, it.next().getHaplotype());
  }
}

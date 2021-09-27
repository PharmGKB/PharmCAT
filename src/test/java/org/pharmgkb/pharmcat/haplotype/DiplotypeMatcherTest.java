package org.pharmgkb.pharmcat.haplotype;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.function.Function;
import java.util.stream.Collectors;
import com.google.common.base.Joiner;
import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import org.pharmgkb.common.util.NoDuplicateMergeFunction;
import org.pharmgkb.pharmcat.definition.model.NamedAllele;
import org.pharmgkb.pharmcat.definition.model.VariantLocus;
import org.pharmgkb.pharmcat.haplotype.model.DiplotypeMatch;
import org.pharmgkb.pharmcat.haplotype.model.HaplotypeMatch;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.core.IsEqual.equalTo;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.fail;


/**
 * JUnit test for {@link DiplotypeMatcher}.
 *
 * @author Mark Woon
 */
class DiplotypeMatcherTest {
  private static VariantLocus[] s_positions;
  private static final SortedSet<NamedAllele> s_haplotypes = new TreeSet<>();


  @BeforeAll
  static void beforeClass() {

    // initialize test variants
    VariantLocus var1 = new VariantLocus("chr1", 1, "g.1T>A");
    VariantLocus var2 = new VariantLocus("chr1", 2, "g.2T>A");
    VariantLocus var3 = new VariantLocus("chr1", 3, "g.3T>A");

    s_positions = new VariantLocus[] { var1, var2, var3 };

    // initialize test haplotypes
    String[] alleles = new String[] { "A", "C", "C" };
    NamedAllele hap1 = new NamedAllele("*1", "*1", alleles, alleles, true);
    hap1.initialize(s_positions);

    alleles = new String[] { "G", null, null };
    NamedAllele hap2 = new NamedAllele("*4a", "*4a", alleles, alleles, false);
    hap2.initialize(s_positions);

    alleles = new String[] { "G", "T", "T" };
    NamedAllele hap3 = new NamedAllele("*4b", "*4b", alleles, alleles, false);
    hap3.initialize(s_positions);

    alleles = new String[] { null, "T", "T" };
    NamedAllele hap4 = new NamedAllele("*17", "*17", alleles, alleles, false);
    hap4.initialize(s_positions);

    s_haplotypes.add(hap1);
    s_haplotypes.add(hap2);
    s_haplotypes.add(hap3);
    s_haplotypes.add(hap4);

    /*
            | 1 | 2 | 3 |
        *1  | A | C | C |
        *4a | G |   |   |
        *4b | G | T | T |
        *17 |   | T | T |
    */
  }


  @Test
  void test1() {

    List<SampleAllele> alleles = Arrays.asList(
        new SampleAllele("chr1", 1, "A", "G", false, Lists.newArrayList("A", "G")),
        new SampleAllele("chr1", 2, "C", "C", false, Lists.newArrayList("C", "T")),
        new SampleAllele("chr1", 3, "C", "C", false, Lists.newArrayList("C", "T"))
    );

    List<DiplotypeMatch> pairMatches = computeHaplotypes(alleles);
    List<String> expectedMatches = Lists.newArrayList("*1/*4a");
    assertDiplotypePairs(expectedMatches, pairMatches);
  }


  @Test
  void test2() {

    List<SampleAllele> alleles = Arrays.asList(
        new SampleAllele("chr1", 1, "A", "G", false, Lists.newArrayList("A", "G")),
        new SampleAllele("chr1", 2, "C", "T", false, Lists.newArrayList("C", "T")),
        new SampleAllele("chr1", 3, "C", "T", false, Lists.newArrayList("C", "T"))
    );

    List<DiplotypeMatch> pairMatches = computeHaplotypes(alleles);
    List<String> expectedMatches = Lists.newArrayList("*1/*4b", "*1/*17", "*1/*4a", "*4a/*17");
    assertDiplotypePairs(expectedMatches, pairMatches);
  }


  @Test
  void test3() {

    List<SampleAllele> alleles = Arrays.asList(
        new SampleAllele("chr1", 1, "A", "A", false, Lists.newArrayList("A", "G")),
        new SampleAllele("chr1", 2, "C", "T", false, Lists.newArrayList("C", "T")),
        new SampleAllele("chr1", 3, "C", "T", false, Lists.newArrayList("C", "T"))
    );

    List<DiplotypeMatch> pairMatches = computeHaplotypes(alleles);
    List<String> expectedMatches = Lists.newArrayList("*1/*17");
    assertDiplotypePairs(expectedMatches, pairMatches);
  }


  @Test
  void test4() {

    List<SampleAllele> alleles = Arrays.asList(
        new SampleAllele("chr1", 1, "G", "G", false, Lists.newArrayList("A", "G")),
        new SampleAllele("chr1", 2, "T", "T", false, Lists.newArrayList("C", "T")),
        new SampleAllele("chr1", 3, "C", "T", false, Lists.newArrayList("C", "T"))
    );

    List<DiplotypeMatch> pairMatches = computeHaplotypes(alleles);
    List<String> expectedMatches = Lists.newArrayList("*4a/*4b", "*4a/*17", "*4a/*4a");
    assertDiplotypePairs(expectedMatches, pairMatches);
  }


  @Test
  void test5() {

    List<SampleAllele> alleles = Arrays.asList(
        new SampleAllele("chr1", 1, "G", "G", false, Lists.newArrayList("A", "G")),
        new SampleAllele("chr1", 2, "C", "T", false, Lists.newArrayList("C", "T")),
        new SampleAllele("chr1", 3, "C", "T", false, Lists.newArrayList("C", "T"))
    );
    List<DiplotypeMatch> pairMatches = computeHaplotypes(alleles);
    List<String> expectedMatches = Lists.newArrayList("*4a/*4b", "*4a/*17", "*4a/*4a");
    assertDiplotypePairs(expectedMatches, pairMatches);
  }


  private List<DiplotypeMatch> computeHaplotypes(List<SampleAllele> alleles) {

    SortedMap<String, SampleAllele> sampleAlleleMap = alleles.stream()
        .collect(Collectors.toMap(s -> "chr1:" + s.getPosition(),
        Function.identity(), new NoDuplicateMergeFunction<>(), TreeMap::new));

    MatchData dataset = new MatchData(sampleAlleleMap, s_positions, null, null);
    dataset.marshallHaplotypes(s_haplotypes);
    dataset.generateSamplePermutations();

    return new DiplotypeMatcher(dataset).compute();
  }



  private void assertDiplotypePairs(List<String> expectedPairs, List<DiplotypeMatch> matches) {

    Preconditions.checkNotNull(expectedPairs);
    Preconditions.checkNotNull(matches);

    List<String> pairs = matches.stream()
        .map(DiplotypeMatch::getName)
        .collect(Collectors.toList());
    assertEquals(matches.size(), new HashSet<>(pairs).size(), "Incoming matches has non-unique pairs");

    if (expectedPairs.size() != pairs.size() || !expectedPairs.equals(pairs)) {
      System.out.println("Expected: [" + Joiner.on(", ").join(expectedPairs));
      System.out.println("Got:      " + pairs);
      fail("Did not get expected matches");
    }
  }


  @Test
  void testComparePermutations() {

    VariantLocus var1 = new VariantLocus("chr1", 1, "g.1T>A");
    VariantLocus var2 = new VariantLocus("chr1", 2, "g.2T>A");
    VariantLocus var3 = new VariantLocus("chr1", 3, "g.3T>A");
    VariantLocus var4 = new VariantLocus("chr1", 4, "g.3T>A");
    VariantLocus[] variants = new VariantLocus[] { var1, var2, var3, var4 };

    String[] alleles = new String[] { "T", "A", "C", "C" };
    NamedAllele hap1 = new NamedAllele("*1", "*1", alleles, alleles, true);
    hap1.initialize(variants);

    alleles = new String[] { null, "T", "C", null };
    NamedAllele hap2 = new NamedAllele("*2", "*2", alleles, alleles, false);
    hap2.initialize(variants);

    alleles = new String[] { null, null, "GG", null };
    NamedAllele hap3 = new NamedAllele("*3", "*3", alleles, alleles, false);
    hap3.initialize(variants);

    Set<String> permutations = Sets.newHashSet(
        "1:T;2:A;3:C;4:C;",
        "1:T;2:A;3:C;4:G;",
        "1:T;2:T;3:C;4:C;",
        "1:T;2:T;3:C;4:G;"
    );
    SortedMap<String, SampleAllele> sampleAlleleMap = new TreeMap<>();
    sampleAlleleMap.put("chr1:1", new SampleAllele("chr1", 1, "T", "T", true, Lists.newArrayList("T")));
    sampleAlleleMap.put("chr1:2", new SampleAllele("chr1", 2, "A", "T", false, Lists.newArrayList("T")));
    sampleAlleleMap.put("chr1:3", new SampleAllele("chr1", 3, "C", "C", false, Lists.newArrayList("C")));
    sampleAlleleMap.put("chr1:4", new SampleAllele("chr1", 4, "C", "G", false, Lists.newArrayList("C")));

    MatchData dataset = new MatchData(sampleAlleleMap, variants, null, null);
    dataset.marshallHaplotypes(new TreeSet<>(Lists.newArrayList(hap1, hap2, hap3)));
    dataset.generateSamplePermutations();
    assertThat(dataset.getPermutations(), equalTo(permutations));

    DiplotypeMatcher diplotypeMatcher = new DiplotypeMatcher(dataset);

    SortedSet<HaplotypeMatch> matches = diplotypeMatcher.comparePermutations();
    assertEquals(2, matches.size());
    Iterator<HaplotypeMatch> it = matches.iterator();
    assertEquals(hap1, it.next().getHaplotype());
    assertEquals(hap2, it.next().getHaplotype());
  }
}

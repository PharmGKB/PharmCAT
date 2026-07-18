package org.pharmgkb.pharmcat.haplotype;

import java.io.IOException;
import java.lang.reflect.Method;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
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
import org.pharmgkb.pharmcat.Env;
import org.pharmgkb.pharmcat.ReportableException;
import org.pharmgkb.pharmcat.definition.model.NamedAllele;
import org.pharmgkb.pharmcat.definition.model.VariantLocus;
import org.pharmgkb.pharmcat.haplotype.model.DiplotypeMatch;
import org.pharmgkb.pharmcat.haplotype.model.HaplotypeMatch;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.core.IsEqual.equalTo;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.fail;


/**
 * JUnit test for {@link DiplotypeMatcher}.
 *
 * @author Mark Woon
 */
class DiplotypeMatcherTest {
  private static Env s_env = null;
  private static VariantLocus[] s_positions;
  private static final SortedSet<NamedAllele> s_haplotypes = new TreeSet<>();


  @BeforeAll
  static void beforeClass() throws IOException, ReportableException {

    s_env = new Env();
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

    SortedSet<SampleAllele> alleles = new TreeSet<>(Arrays.asList(
        new SampleAllele("chr1", 1, "A", "G", false, Lists.newArrayList("A", "G"), "0/1"),
        new SampleAllele("chr1", 2, "C", "C", false, Lists.newArrayList("C", "T"), "0/0"),
        new SampleAllele("chr1", 3, "C", "C", false, Lists.newArrayList("C", "T"), "0/0")
    ));

    SortedSet<DiplotypeMatch> pairMatches = computeHaplotypes(alleles);
    List<String> expectedMatches = Lists.newArrayList("*1/*4a");
    assertDiplotypePairs(expectedMatches, pairMatches);
  }


  @Test
  void test2() {

    SortedSet<SampleAllele> alleles = new TreeSet<>(Arrays.asList(
        new SampleAllele("chr1", 1, "A", "G", false, Lists.newArrayList("A", "G"), "0/1"),
        new SampleAllele("chr1", 2, "C", "T", false, Lists.newArrayList("C", "T"), "0/1"),
        new SampleAllele("chr1", 3, "C", "T", false, Lists.newArrayList("C", "T"), "0/1")
    ));

    SortedSet<DiplotypeMatch> pairMatches = computeHaplotypes(alleles);
    List<String> expectedMatches = Lists.newArrayList("*1/*4b", "*1/*17", "*1/*4a", "*4a/*17");
    assertDiplotypePairs(expectedMatches, pairMatches);
  }


  @Test
  void test3() {

    SortedSet<SampleAllele> alleles = new TreeSet<>(Arrays.asList(
        new SampleAllele("chr1", 1, "A", "A", false, Lists.newArrayList("A", "G"), "0/0"),
        new SampleAllele("chr1", 2, "C", "T", false, Lists.newArrayList("C", "T"), "0/1"),
        new SampleAllele("chr1", 3, "C", "T", false, Lists.newArrayList("C", "T"), "0/1")
    ));

    SortedSet<DiplotypeMatch> pairMatches = computeHaplotypes(alleles);
    List<String> expectedMatches = Lists.newArrayList("*1/*17");
    assertDiplotypePairs(expectedMatches, pairMatches);
  }


  @Test
  void test4() {

    SortedSet<SampleAllele> alleles = new TreeSet<>(Arrays.asList(
        new SampleAllele("chr1", 1, "G", "G", false, Lists.newArrayList("A", "G"), "1/1"),
        new SampleAllele("chr1", 2, "T", "T", false, Lists.newArrayList("C", "T"), "1/1"),
        new SampleAllele("chr1", 3, "C", "T", false, Lists.newArrayList("C", "T"), "0/1")
    ));

    SortedSet<DiplotypeMatch> pairMatches = computeHaplotypes(alleles);
    List<String> expectedMatches = Lists.newArrayList("*4a/*4b", "*4a/*17", "*4a/*4a");
    assertDiplotypePairs(expectedMatches, pairMatches);
  }


  @Test
  void test5() {

    SortedSet<SampleAllele> alleles = new TreeSet<>(Arrays.asList(
        new SampleAllele("chr1", 1, "G", "G", false, Lists.newArrayList("A", "G"), "1/1"),
        new SampleAllele("chr1", 2, "C", "T", false, Lists.newArrayList("C", "T"), "0/1"),
        new SampleAllele("chr1", 3, "C", "T", false, Lists.newArrayList("C", "T"), "0/1")
    ));
    SortedSet<DiplotypeMatch> pairMatches = computeHaplotypes(alleles);
    List<String> expectedMatches = Lists.newArrayList("*4a/*4b", "*4a/*17", "*4a/*4a");
    assertDiplotypePairs(expectedMatches, pairMatches);
  }


  @Test
  void testStructuredComplementCheckWithNullAllele() {

    VariantLocus var1 = new VariantLocus("chr1", 1, "g.1T>A");
    VariantLocus var2 = new VariantLocus("chr1", 2, "g.2T>A");
    var1.setRef("A");
    var2.setRef("C");
    VariantLocus[] variants = new VariantLocus[] { var1, var2 };

    String[] alleles = new String[] { "A", "C" };
    NamedAllele hap1 = new NamedAllele("*1", "*1", alleles, alleles, true);
    hap1.initialize(variants);

    alleles = new String[] { null, "T" };
    NamedAllele hap2 = new NamedAllele("*2", "*2", alleles, alleles, false);
    hap2.initialize(variants);

    SortedMap<String, SampleAllele> sampleAlleleMap = new TreeMap<>();
    sampleAlleleMap.put("chr1:1", new SampleAllele("chr1", 1, "A", null, false,
        Lists.newArrayList("A"), "0"));
    sampleAlleleMap.put("chr1:2", new SampleAllele("chr1", 2, "C", "T", false,
        Lists.newArrayList("C", "T"), "0/1"));

    MatchData dataset = new MatchData("Sample_1", "CYP2B6", sampleAlleleMap, variants, null, null);
    dataset.marshallHaplotypes("TEST", new TreeSet<>(Lists.newArrayList(hap1, hap2)), false);
    dataset.generateSamplePermutations();

    SortedSet<DiplotypeMatch> matches = new DiplotypeMatcher(s_env, dataset)
        .compute(false, false);

    assertEquals(1, matches.size());
    assertEquals("*1/*2", matches.first().getName());
  }


  @Test
  void testLegacyComplementCheckWithMissingSequencePosition() throws Exception {

    VariantLocus var1 = new VariantLocus("chr1", 1, "g.1T>A");
    VariantLocus var2 = new VariantLocus("chr1", 2, "g.2T>A");
    var1.setRef("A");
    var2.setRef("C");
    VariantLocus[] variants = new VariantLocus[] { var1, var2 };

    String[] alleles = new String[] { "A", "C" };
    NamedAllele hap = new NamedAllele("*1", "*1", alleles, alleles, true);
    hap.initialize(variants);

    SortedMap<String, SampleAllele> sampleAlleleMap = new TreeMap<>();
    sampleAlleleMap.put("chr1:1", new SampleAllele("chr1", 1, "A", "A", false,
        Lists.newArrayList("A"), "0/0"));
    sampleAlleleMap.put("chr1:2", new SampleAllele("chr1", 2, "C", "T", false,
        Lists.newArrayList("C", "T"), "0/1"));

    MatchData dataset = new MatchData("Sample_1", "CYP2B6", sampleAlleleMap, variants, null, null);
    dataset.marshallHaplotypes("TEST", new TreeSet<>(Lists.newArrayList(hap)), false);
    dataset.generateSamplePermutations();

    Method method = DiplotypeMatcher.class.getDeclaredMethod("isViableComplement", String.class, String.class);
    method.setAccessible(true);

    boolean viable = (boolean)method.invoke(new DiplotypeMatcher(s_env, dataset), "1:A", "1:A");
    assertFalse(viable);
  }


  @Test
  void testStructuredScoringWithWobble() {

    VariantLocus var1 = new VariantLocus("chr1", 1, "g.1C>T");
    VariantLocus var2 = new VariantLocus("chr1", 2, "g.2A>G");
    var1.setRef("C");
    var2.setRef("A");
    VariantLocus[] variants = new VariantLocus[] { var1, var2 };

    String[] alleles = new String[] { "C", "A" };
    NamedAllele ref = new NamedAllele("*1", "*1", alleles, alleles, true);
    ref.initialize(variants);

    alleles = new String[] { "Y", null };
    NamedAllele wobble = new NamedAllele("*2", "*2", alleles, alleles, false);
    wobble.initialize(variants);

    SortedMap<String, SampleAllele> sampleAlleleMap = new TreeMap<>();
    sampleAlleleMap.put("chr1:1", new SampleAllele("chr1", 1, "C", "C", false,
        Lists.newArrayList("C", "T"), "0/0"));
    sampleAlleleMap.put("chr1:2", new SampleAllele("chr1", 2, "A", "A", false,
        Lists.newArrayList("A", "G"), "0/0"));

    MatchData dataset = new MatchData("Sample_1", "CYP2B6", sampleAlleleMap, variants, null, null);
    dataset.marshallHaplotypes("TEST", new TreeSet<>(Lists.newArrayList(ref, wobble)), false);
    dataset.generateSamplePermutations();

    SortedSet<DiplotypeMatch> matches = new DiplotypeMatcher(s_env, dataset)
        .compute(false, false);
    Map<String, Integer> scores = matches.stream()
        .collect(Collectors.toMap(DiplotypeMatch::getName, DiplotypeMatch::getScore));
    assertEquals(4, scores.get("*1/*1"));
    assertEquals(2, scores.get("*1/*2"));
    assertEquals(0, scores.get("*2/*2"));

    SortedSet<DiplotypeMatch> topMatches = new DiplotypeMatcher(s_env, dataset)
        .compute(false, true);
    assertEquals(1, topMatches.size());
    assertEquals("*1/*1", topMatches.first().getName());
  }


  private SortedSet<DiplotypeMatch> computeHaplotypes(SortedSet<SampleAllele> alleles) {

    SortedMap<String, SampleAllele> sampleAlleleMap = alleles.stream()
        .collect(Collectors.toMap(s -> "chr1:" + s.getPosition(),
        Function.identity(), new NoDuplicateMergeFunction<>(), TreeMap::new));

    MatchData dataset = new MatchData("Sample_1", "CYP2B6", sampleAlleleMap, s_positions, null, null);
    dataset.marshallHaplotypes("TEST", s_haplotypes, false);
    dataset.generateSamplePermutations();

    return new DiplotypeMatcher(s_env, dataset)
        .compute(false, false);
  }



  private void assertDiplotypePairs(List<String> expectedPairs, SortedSet<DiplotypeMatch> matches) {

    Preconditions.checkNotNull(expectedPairs);
    Preconditions.checkNotNull(matches);

    List<String> pairs = matches.stream()
        .map(DiplotypeMatch::getName)
        .toList();
    assertEquals(matches.size(), new HashSet<>(pairs).size(), "Incoming matches has non-unique pairs");

    if (expectedPairs.size() != pairs.size() || !expectedPairs.equals(pairs)) {
      System.out.println("Expected: [" + Joiner.on(", ").join(expectedPairs));
      System.out.println("Got:      " + pairs);
      fail("Did not get expected matches");
    }
  }


  @Test
  void testMarshallHaplotypesDropsDefinitionsWithOnlyMissingPositions() {

    VariantLocus var1 = new VariantLocus("chr1", 1, "g.1A>G");
    VariantLocus var2 = new VariantLocus("chr1", 2, "g.2C>T");
    VariantLocus var3 = new VariantLocus("chr1", 3, "g.3C>T");
    VariantLocus[] variants = new VariantLocus[] { var1, var2, var3 };

    NamedAllele ref = new NamedAllele("*1", "*1", new String[] { "A", "C", "C" },
        new String[] { "A", "C", "C" }, true);
    ref.initialize(variants);
    NamedAllele available = new NamedAllele("*2", "*2", new String[] { "G", null, null },
        new String[] { "G", null, null }, false);
    available.initialize(variants);
    NamedAllele missing = new NamedAllele("*3", "*3", new String[] { null, "T", "T" },
        new String[] { null, "T", "T" }, false);
    missing.initialize(variants);

    SortedMap<String, SampleAllele> sampleAlleleMap = new TreeMap<>();
    sampleAlleleMap.put("chr1:1", new SampleAllele("chr1", 1, "A", "G", false,
        Lists.newArrayList("A", "G"), "0/1"));

    MatchData dataset = new MatchData("Sample_1", "GENE", sampleAlleleMap, variants, null, null);
    dataset.marshallHaplotypes("GENE", new TreeSet<>(List.of(ref, available, missing)), false);

    assertEquals(List.of("*1", "*2"), dataset.getHaplotypes().stream()
        .map(NamedAllele::getName)
        .toList());
  }


  @Test
  void testComparePermutations() {

    VariantLocus var1 = new VariantLocus("chr1", 1, "g.1T>A");
    VariantLocus var2 = new VariantLocus("chr1", 2, "g.2T>A");
    VariantLocus var3 = new VariantLocus("chr1", 3, "g.3T>A");
    VariantLocus var4 = new VariantLocus("chr1", 4, "g.3T>A");
    // Keep definition positions out of order to exercise position-sorted internal permutation matching.
    VariantLocus[] variants = new VariantLocus[] { var3, var1, var4, var2 };

    String[] alleles = new String[] { "C", "T", "C", "A" };
    NamedAllele hap1 = new NamedAllele("*1", "*1", alleles, alleles, true);
    hap1.initialize(variants);

    alleles = new String[] { "C", null, null, "T" };
    NamedAllele hap2 = new NamedAllele("*2", "*2", alleles, alleles, false);
    hap2.initialize(variants);

    alleles = new String[] { "GG", null, null, null };
    NamedAllele hap3 = new NamedAllele("*3", "*3", alleles, alleles, false);
    hap3.initialize(variants);

    Set<String> permutations = Sets.newHashSet(
        "1:T;2:A;3:C;4:C",
        "1:T;2:A;3:C;4:G",
        "1:T;2:T;3:C;4:C",
        "1:T;2:T;3:C;4:G"
    );
    SortedMap<String, SampleAllele> sampleAlleleMap = new TreeMap<>();
    sampleAlleleMap.put("chr1:1", new SampleAllele("chr1", 1, "T", "T", true, Lists.newArrayList("T"), "0/0"));
    sampleAlleleMap.put("chr1:2", new SampleAllele("chr1", 2, "A", "T", false, Lists.newArrayList("T"), "1/0"));
    sampleAlleleMap.put("chr1:3", new SampleAllele("chr1", 3, "C", "C", false, Lists.newArrayList("C"), "0/0"));
    sampleAlleleMap.put("chr1:4", new SampleAllele("chr1", 4, "C", "G", false, Lists.newArrayList("C"), "0/1"));

    MatchData dataset = new MatchData("Sample_1", "GENE", sampleAlleleMap, variants, null, null);
    dataset.marshallHaplotypes("TEST", new TreeSet<>(Lists.newArrayList(hap1, hap2, hap3)), false);
    dataset.generateSamplePermutations();
    assertThat(dataset.getPermutationStrings(), equalTo(permutations));

    SortedSet<HaplotypeMatch> matches = dataset.comparePermutations();
    assertEquals(2, matches.size());
    for (HaplotypeMatch match : matches) {
      for (String sequence : match.getSequences()) {
        assertNotNull(match.getSequencePermutation(sequence));
      }
    }
    Iterator<HaplotypeMatch> it = matches.iterator();
    assertEquals(hap1, it.next().getHaplotype());
    assertEquals(hap2, it.next().getHaplotype());
  }


  @Test
  void testComparePermutationsWithAmbiguousAndRepeatAlleles() {

    VariantLocus var1 = new VariantLocus("chr1", 1, "g.1T>A");
    VariantLocus var2 = new VariantLocus("chr1", 2, "g.2T>A");
    VariantLocus[] variants = new VariantLocus[] { var1, var2 };

    NamedAllele ambiguous = new NamedAllele("*1", "*1", new String[] { "Y", null },
        new String[] { "Y", null }, false);
    ambiguous.initialize(variants);
    NamedAllele repeat = new NamedAllele("*2", "*2", new String[] { null, "CAT or CATAT" },
        new String[] { null, "CAT or CATAT" }, false);
    repeat.initialize(variants);

    SortedMap<String, SampleAllele> sampleAlleleMap = new TreeMap<>();
    sampleAlleleMap.put("chr1:1", new SampleAllele("chr1", 1, "C", "T", false,
        Lists.newArrayList("C", "T"), "0/1"));
    sampleAlleleMap.put("chr1:2", new SampleAllele("chr1", 2, "CAT", "CATAT", false,
        Lists.newArrayList("CAT", "CATAT"), "0/1"));

    MatchData dataset = new MatchData("Sample_1", "GENE", sampleAlleleMap, variants, null, null);
    dataset.marshallHaplotypes("TEST", new TreeSet<>(List.of(ambiguous, repeat)), false);
    dataset.generateSamplePermutations();

    SortedSet<HaplotypeMatch> matches = dataset.comparePermutations();
    assertEquals(2, matches.size());
    Iterator<HaplotypeMatch> it = matches.iterator();
    assertEquals(ambiguous, it.next().getHaplotype());
    assertEquals(repeat, it.next().getHaplotype());
  }


  @Test
  void testComparePermutationsAfterRemarshallingHaplotypes() {

    VariantLocus variant = new VariantLocus("chr1", 1, "g.1T>A");
    VariantLocus[] variants = new VariantLocus[] { variant };
    NamedAllele cAllele = new NamedAllele("1", "*1", new String[] { "C" }, new String[] { "C" }, false);
    cAllele.initialize(variants);
    NamedAllele tAllele = new NamedAllele("2", "*2", new String[] { "T" }, new String[] { "T" }, false);
    tAllele.initialize(variants);

    SortedMap<String, SampleAllele> sampleAlleleMap = new TreeMap<>();
    sampleAlleleMap.put("chr1:1", new SampleAllele("chr1", 1, "C", "C", false,
        Lists.newArrayList("C", "T"), "0/0"));
    MatchData dataset = new MatchData("Sample_1", "GENE", sampleAlleleMap, variants, null, null);
    dataset.marshallHaplotypes("TEST", new TreeSet<>(List.of(cAllele)), false);
    dataset.generateSamplePermutations();
    assertEquals(cAllele, dataset.comparePermutations().first().getHaplotype());

    dataset.marshallHaplotypes("TEST", new TreeSet<>(List.of(tAllele)), false);
    assertEquals(0, dataset.comparePermutations().size());
  }


  @Test
  void testLazyReferenceDefaultingMaterializesMatchedHaplotype() {

    VariantLocus var1 = new VariantLocus("chr1", 1, "g.1A>G");
    VariantLocus var2 = new VariantLocus("chr1", 2, "g.2C>T");
    var1.setRef("A");
    var2.setRef("C");
    VariantLocus[] variants = new VariantLocus[] { var1, var2 };

    String[] alleles = new String[] { "A", "C" };
    NamedAllele ref = new NamedAllele("*1", "*1", alleles, alleles, true);
    ref.initialize(variants);

    alleles = new String[] { "G", null };
    NamedAllele hap = new NamedAllele("*2", "*2", alleles, alleles, false);
    hap.initialize(variants);

    SortedMap<String, SampleAllele> sampleAlleleMap = new TreeMap<>();
    sampleAlleleMap.put("chr1:1", new SampleAllele("chr1", 1, "G", "G", false,
        Lists.newArrayList("A", "G"), "1/1"));
    sampleAlleleMap.put("chr1:2", new SampleAllele("chr1", 2, "C", "C", false,
        Lists.newArrayList("C", "T"), "0/0"));

    MatchData dataset = new MatchData("Sample_1", "GENE", sampleAlleleMap, variants, null, null);
    dataset.marshallHaplotypes("TEST", new TreeSet<>(Lists.newArrayList(ref, hap)), false);
    dataset.defaultMissingAllelesToReference();
    dataset.generateSamplePermutations();

    SortedSet<HaplotypeMatch> matches = dataset.comparePermutations();
    assertEquals(1, matches.size());
    NamedAllele matchedHaplotype = matches.first().getHaplotype();
    assertEquals("*2", matchedHaplotype.getName());
    assertEquals("G", matchedHaplotype.getAllele(var1));
    assertEquals("C", matchedHaplotype.getAllele(var2));

    NamedAllele internalHaplotype = dataset.getHaplotypes().stream()
        .filter(h -> h.getName().equals("*2"))
        .findAny()
        .orElseThrow();
    assertNull(internalHaplotype.getAllele(var2));
    NamedAllele outputHaplotype = dataset.getHaplotypesForOutput().stream()
        .filter(h -> h.getName().equals("*2"))
        .findAny()
        .orElseThrow();
    assertEquals("C", outputHaplotype.getAllele(var2));
  }


  @Test
  void testCombinationSimpleMatchRetainsStructuredPermutation() {

    VariantLocus var1 = new VariantLocus("chr1", 1, "g.1T>A");
    VariantLocus var2 = new VariantLocus("chr1", 2, "g.2T>A");
    var1.setRef("A");
    var2.setRef("C");
    VariantLocus[] variants = new VariantLocus[] { var1, var2 };

    String[] alleles = new String[] { "A", "C" };
    NamedAllele ref = new NamedAllele("*1", "*1", alleles, alleles, true);
    ref.initialize(variants);

    alleles = new String[] { "G", null };
    NamedAllele hap = new NamedAllele("*2", "*2", alleles, alleles, false);
    hap.initialize(variants);

    SortedMap<String, SampleAllele> sampleAlleleMap = new TreeMap<>();
    sampleAlleleMap.put("chr1:1", new SampleAllele("chr1", 1, "G", "G", false,
        Lists.newArrayList("A", "G"), "1/1"));
    sampleAlleleMap.put("chr1:2", new SampleAllele("chr1", 2, "C", "C", false,
        Lists.newArrayList("C"), "0/0"));

    MatchData dataset = new MatchData("Sample_1", "CYP2B6", sampleAlleleMap, variants, null, null);
    dataset.marshallHaplotypes("TEST", new TreeSet<>(Lists.newArrayList(ref, hap)), true);
    dataset.generateSamplePermutations();

    SortedSet<HaplotypeMatch> matches = new CombinationMatcher(s_env.getDefinitionReader().getDefinitionFile("CYP2B6"),
        false).compute(dataset).stream()
        .map(m -> (HaplotypeMatch)m)
        .collect(Collectors.toCollection(TreeSet::new));

    assertEquals(1, matches.size());
    HaplotypeMatch match = matches.first();
    assertEquals("*2", match.getName());
    for (String sequence : match.getSequences()) {
      assertNotNull(match.getSequencePermutation(sequence));
    }
  }


  @Test
  void testCombinationCandidateIndexWithWobbleCorePosition() {

    VariantLocus var1 = new VariantLocus("chr1", 1, "g.1C>T");
    VariantLocus var2 = new VariantLocus("chr1", 2, "g.2A>G");
    var1.setRef("C");
    var2.setRef("A");
    VariantLocus[] variants = new VariantLocus[] { var1, var2 };

    String[] alleles = new String[] { "C", "A" };
    NamedAllele ref = new NamedAllele("*1", "*1", alleles, alleles, true);
    ref.initialize(variants);

    alleles = new String[] { "Y", null };
    NamedAllele wobble = new NamedAllele("*2", "*2", alleles, alleles, false);
    wobble.initialize(variants);

    SortedMap<String, SampleAllele> sampleAlleleMap = new TreeMap<>();
    sampleAlleleMap.put("chr1:1", new SampleAllele("chr1", 1, "T", "T", false,
        Lists.newArrayList("C", "T"), "1/1"));
    sampleAlleleMap.put("chr1:2", new SampleAllele("chr1", 2, "A", "A", false,
        Lists.newArrayList("A", "G"), "0/0"));

    MatchData dataset = new MatchData("Sample_1", "CYP2B6", sampleAlleleMap, variants, null, null);
    dataset.marshallHaplotypes("TEST", new TreeSet<>(Lists.newArrayList(ref, wobble)), true);
    dataset.generateSamplePermutations();

    SortedSet<HaplotypeMatch> matches = new CombinationMatcher(s_env.getDefinitionReader().getDefinitionFile("CYP2B6"),
        false).compute(dataset).stream()
        .map(m -> (HaplotypeMatch)m)
        .collect(Collectors.toCollection(TreeSet::new));

    assertEquals(1, matches.size());
    assertEquals("*2", matches.first().getName());
  }
}

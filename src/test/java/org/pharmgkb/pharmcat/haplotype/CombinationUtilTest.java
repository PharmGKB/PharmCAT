package org.pharmgkb.pharmcat.haplotype;

import java.util.Arrays;
import java.util.List;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import org.junit.jupiter.api.Test;
import org.pharmgkb.pharmcat.definition.model.NamedAllele;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;


/**
 * JUnit test for {@link CombinationUtil}.
 *
 * @author Mark Woon
 */
class CombinationUtilTest {


  @Test
  void testGeneratePermutationsNotPhased() {

    List<SampleAllele> alleles = Arrays.asList(
        new SampleAllele("chr1", 1, "T", "T", false, Lists.newArrayList("T", "C"), "0/0"),
        new SampleAllele("chr1", 2, "A", "T", false, Lists.newArrayList("A", "T"),  "0/1"),
        new SampleAllele("chr1", 3, "C", "C", false, Lists.newArrayList("C", "C"), "0/0"),
        new SampleAllele("chr1", 4, "C", "G", false, Lists.newArrayList("C", "G"), "0/1")
    );

    Set<String> expectedPermutations = Sets.newHashSet(
        "1:T;2:A;3:C;4:C;",
        "1:T;2:A;3:C;4:G;",
        "1:T;2:T;3:C;4:C;",
        "1:T;2:T;3:C;4:G;"
        );
    Set<String> permutations = CombinationUtil.generatePermutations(alleles);
    assertEquals(expectedPermutations.size(), permutations.size());
    for (String p : permutations) {
      assertTrue(expectedPermutations.remove(p));
    }
  }

  @Test
  void testGeneratePermutationPhased() {

    List<SampleAllele> alleles = Arrays.asList(
        new SampleAllele("chr1", 1, "T", "T", true, Lists.newArrayList("T", "C"), "0|0"),
        new SampleAllele("chr1", 2, "A", "T", true, Lists.newArrayList("A", "T"), "0|1"),
        new SampleAllele("chr1", 3, "C", "C", true, Lists.newArrayList("C", "C"), "0|0"),
        new SampleAllele("chr1", 4, "C", "G", true, Lists.newArrayList("C", "G"), "0|1")
    );

    Set<String> expectedPermutations = Sets.newHashSet(
        "1:T;2:A;3:C;4:C;",
        "1:T;2:T;3:C;4:G;"
    );
    Set<String> permutations = CombinationUtil.generatePermutations(alleles);
    assertEquals(expectedPermutations.size(), permutations.size());
    for (String p : permutations) {
      assertTrue(expectedPermutations.remove(p));
    }
  }


  @Test
  void testGeneratePermutationPhaseSets() {

    List<SampleAllele> alleles = Arrays.asList(
        new SampleAllele("chr1", 1, "T", "T", true, true, null, Lists.newArrayList("T", "C"), "0|0", null, false),
        new SampleAllele("chr1", 2, "A", "T", true, false, null, Lists.newArrayList("A", "T"), "0|1", null, false),
        new SampleAllele("chr1", 3, "C", "C", true, true, 1, Lists.newArrayList("C", "C"), "0|0", null, false),
        new SampleAllele("chr1", 4, "C", "G", true, false, 1, Lists.newArrayList("C", "G"), "0|1", null, false)
    );

    Set<String> expectedPermutations = Sets.newHashSet(
        "1:T;2:A;3:C;4:C;",
        "1:T;2:T;3:C;4:G;",
        "1:T;2:A;3:C;4:G;",
        "1:T;2:T;3:C;4:C;"
    );
    Set<String> permutations = CombinationUtil.generatePermutations(alleles);
    System.out.println(permutations);
    assertEquals(expectedPermutations.size(), permutations.size());
    for (String p : permutations) {
      assertTrue(expectedPermutations.remove(p));
    }
  }


  @Test
  void testGeneratePermutationPhaseSetsWithUnphased() {

    List<SampleAllele> alleles = Arrays.asList(
        new SampleAllele("chr1", 1, "T", "C", true, false, 2, Lists.newArrayList("T", "C"), "0|1", null, false),
        new SampleAllele("chr1", 2, "A", "T", true, false, 2, Lists.newArrayList("A", "T"), "0|1", null, false),
        new SampleAllele("chr1", 3, "G", "C", false, false, 1, Lists.newArrayList("G", "C"), "1/0", null, false),
        new SampleAllele("chr1", 4, "C", "C", true, true, 1, Lists.newArrayList("C", "C"), "0|0", null, false),
        new SampleAllele("chr1", 5, "C", "G", true, false, 1, Lists.newArrayList("C", "G"), "0|1", null, false)
    );

    Set<String> expectedPermutations = Sets.newHashSet(
        "1:T;2:A;3:G;4:C;5:C;",
        "1:T;2:A;3:C;4:C;5:C;",

        "1:T;2:A;3:G;4:C;5:G;",
        "1:T;2:A;3:C;4:C;5:G;",

        "1:C;2:T;3:G;4:C;5:C;",
        "1:C;2:T;3:C;4:C;5:C;",

        "1:C;2:T;3:G;4:C;5:G;",
        "1:C;2:T;3:C;4:C;5:G;"

        );
    Set<String> permutations = CombinationUtil.generatePermutations(alleles);
    permutations.forEach(System.out::println);
    assertEquals(expectedPermutations.size(), permutations.size());
    for (String p : permutations) {
      assertTrue(expectedPermutations.remove(p));
    }
  }


  @Test
  void testGeneratePerfectPairs() {

    SortedSet<NamedAllele> haplotypes = new TreeSet<>();
    String[] alleles = new String[] { "T" };

    NamedAllele hap = new NamedAllele("*1", "*1", alleles, alleles, true);
    haplotypes.add(hap);

    hap = new NamedAllele("*4", "*4", alleles, alleles, false);
    haplotypes.add(hap);

    hap = new NamedAllele("*3", "*3", alleles, alleles, false);
    haplotypes.add(hap);

    hap = new NamedAllele("*2", "*2", alleles, alleles, false);
    haplotypes.add(hap);

    assertEquals(4, haplotypes.size());

    Set<String> expectedPairs = Sets.newHashSet(
        "*1*1",
        "*1*2",
        "*1*3",
        "*1*4",
        "*2*2",
        "*2*3",
        "*2*4",
        "*3*3",
        "*3*4",
        "*4*4"
    );
    List<List<NamedAllele>> pairs = CombinationUtil.generatePerfectPairs(haplotypes);
    assertEquals(expectedPairs.size(), pairs.size());

    for (List<NamedAllele> pair : pairs) {
      String p = pair.get(0).getName() + pair.get(1).getName();
      assertTrue(expectedPairs.remove(p));
    }
    assertEquals(0, expectedPairs.size());
  }
}

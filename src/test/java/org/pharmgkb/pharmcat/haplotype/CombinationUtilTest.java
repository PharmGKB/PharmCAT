package org.pharmgkb.pharmcat.haplotype;

import java.util.Arrays;
import java.util.List;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import org.junit.Test;
import org.pharmgkb.pharmcat.definition.model.NamedAllele;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;


/**
 * JUnit test for {@link CombinationUtil}.
 *
 * @author Mark Woon
 */
public class CombinationUtilTest {


  @Test
  public void testGeneratePermutationsNotPhased() {

    List<SampleAllele> alleles = Arrays.asList(
        new SampleAllele("chr1", 1, "T", "T", false, Lists.newArrayList("T", "C")),
        new SampleAllele("chr1", 2, "A", "T", false, Lists.newArrayList("A", "T")),
        new SampleAllele("chr1", 3, "C", "C", false, Lists.newArrayList("C", "C")),
        new SampleAllele("chr1", 4, "C", "G", false, Lists.newArrayList("C", "G"))
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
      System.out.println(p);
      assertTrue(expectedPermutations.remove(p));
    }
  }

  @Test
  public void testGeneratePermutationPhased() {

    List<SampleAllele> alleles = Arrays.asList(
        new SampleAllele("chr1", 1, "T", "T", true, Lists.newArrayList("T", "C")),
        new SampleAllele("chr1", 2, "A", "T", true, Lists.newArrayList("A", "T")),
        new SampleAllele("chr1", 3, "C", "C", true, Lists.newArrayList("C", "C")),
        new SampleAllele("chr1", 4, "C", "G", true, Lists.newArrayList("C", "G"))
    );

    Set<String> expectedPermutations = Sets.newHashSet(
        "1:T;2:A;3:C;4:C;",
        "1:T;2:T;3:C;4:G;"
    );
    Set<String> permutations = CombinationUtil.generatePermutations(alleles);
    assertEquals(expectedPermutations.size(), permutations.size());
    for (String p : permutations) {
      System.out.println(p);
      assertTrue(expectedPermutations.remove(p));
    }
  }


  @Test
  public void testGeneratePerfectPairs() {

    SortedSet<NamedAllele> haplotypes = new TreeSet<>();

    NamedAllele hap = new NamedAllele("*1", "*1", new String[1]);
    haplotypes.add(hap);

    hap = new NamedAllele("*4", "*4", new String[1]);
    haplotypes.add(hap);

    hap = new NamedAllele("*3", "*3", new String[1]);
    haplotypes.add(hap);

    hap = new NamedAllele("*2", "*2", new String[1]);
    haplotypes.add(hap);

    assertEquals(4, haplotypes.size());
    System.out.println(haplotypes);

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
      System.out.println(pair);
      String p = pair.get(0).getName() + pair.get(1).getName();
      assertTrue(expectedPairs.remove(p));
    }
    assertEquals(0, expectedPairs.size());
  }
}

package org.cpic.haplotype;

import java.util.Arrays;
import java.util.List;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import com.google.common.collect.Sets;
import org.junit.Test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;


/**
 * @author Mark Woon
 */
public class CombinationUtilTest {


  @Test
  public void testGeneratePermutations() {

    List<SampleAllele> alleles = Arrays.asList(
        new SampleAllele("chr1", 1, "T", "T", false),
        new SampleAllele("chr1", 2, "A", "T", false),
        new SampleAllele("chr1", 3, "C", "C", false),
        new SampleAllele("chr1", 4, "C", "G", false)
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
  public void testGeneratePerfectPairs() {

    SortedSet<Haplotype> haplotypes = new TreeSet<>();

    Haplotype hap = new Haplotype();
    hap.setCommonName("*1");
    haplotypes.add(hap);

    hap = new Haplotype();
    hap.setCommonName("*4");
    haplotypes.add(hap);

    hap = new Haplotype();
    hap.setCommonName("*3");
    haplotypes.add(hap);

    hap = new Haplotype();
    hap.setCommonName("*2");
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
    List<List<Haplotype>> pairs = CombinationUtil.generatePerfectPairs(haplotypes);
    assertEquals(expectedPairs.size(), pairs.size());

    for (List<Haplotype> pair : pairs) {
      System.out.println(pair);
      String p = pair.get(0).getCommonName() + pair.get(1).getCommonName();
      assertTrue(expectedPairs.remove(p));
    }
    assertEquals(0, expectedPairs.size());
  }
}

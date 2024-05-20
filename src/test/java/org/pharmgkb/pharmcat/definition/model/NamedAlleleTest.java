package org.pharmgkb.pharmcat.definition.model;

import java.util.Arrays;
import java.util.Collections;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.regex.Pattern;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;


/**
 * JUnit test for {@link NamedAllele}.
 *
 * @author Mark Woon
 */
class NamedAlleleTest {

  @Test
  void testCalculatePermutations() {

    VariantLocus var1 = new VariantLocus("chr1", 1, "g.1T>A");
    VariantLocus var2 = new VariantLocus("chr1", 2, "g.2T>A");
    VariantLocus var3 = new VariantLocus("chr1", 3, "g.3T>A");
    VariantLocus var4 = new VariantLocus("chr1", 4, "g.3T>A");
    VariantLocus[] variants = new VariantLocus[] { var1, var2, var3, var4 };


    String[] alleles = new String[] { "T", "A", "G", "C" };
    NamedAllele ref = new NamedAllele("*1", "*1", alleles, alleles, true, false);
    ref.initialize(variants);

    alleles = new String[] { "C", "T", null, null };
    NamedAllele hap2 = new NamedAllele("*2", "*2", alleles, alleles, false, false);
    hap2.initialize(variants);

    alleles = new String[] { "C", null, null, "Y" };
    NamedAllele hap3 = new NamedAllele("*3", "*3", alleles, alleles, false, false);
    hap3.initialize(variants);

    // permutations should be consistent even if positions are out-of-order
    alleles = new String[] { "Y", null, null, "C" };
    NamedAllele hap4 = new NamedAllele("*4", "*4", alleles, alleles, false, false);
    Arrays.sort(variants, Collections.reverseOrder());
    hap4.initialize(variants);

    assertEquals("1:T;2:A;3:G;4:C;", ref.getPermutations().pattern());
    assertEquals("1:C;2:T;3:.*?;4:.*?;", hap2.getPermutations().pattern());
    assertEquals("1:C;2:.*?;3:.*?;4:[CT];", hap3.getPermutations().pattern());
    assertEquals("1:C;2:.*?;3:.*?;4:[CT];", hap4.getPermutations().pattern());
  }


  @Test
  void testPermutationPattern() {

    String seq = "1:C;2:C;3:C;4:C;";
    Pattern p = Pattern.compile("1:C;2:.?;3:.?;4:[CT];");
    assertTrue(p.matcher(seq).matches());
  }


  @Test
  void testSorting() {

    String[] alleles = new String[] {"A"};
    NamedAllele nonRef = new NamedAllele("1", "*1", alleles, alleles, false, false);
    NamedAllele ref    = new NamedAllele("2", "*2", alleles, alleles, true, false);
    SortedSet<NamedAllele> set = new TreeSet<>();
    set.add(nonRef);
    set.add(ref);
    assertEquals(ref, set.first());
  }
}

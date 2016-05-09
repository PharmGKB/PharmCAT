package org.pharmgkb.pharmcat.definition.model;

import java.util.Arrays;
import java.util.Collections;
import java.util.regex.Pattern;
import org.junit.Test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;


/**
 * JUnit test for {@link NamedAllele}.
 *
 * @author Mark Woon
 */
public class NamedAlleleTest {

  @Test
  public void testCalculatePermutations() {

    VariantLocus var1 = new VariantLocus(1, "g.1T>A");
    VariantLocus var2 = new VariantLocus(2, "g.2T>A");
    VariantLocus var3 = new VariantLocus(3, "g.3T>A");
    VariantLocus var4 = new VariantLocus(4, "g.3T>A");
    VariantLocus[] variants = new VariantLocus[] { var1, var2, var3, var4 };


    NamedAllele ref = new NamedAllele("*1", "*1", new String[] { "T", "A", "G", "C" });
    ref.finalize(variants);

    NamedAllele hap2 = new NamedAllele("*2", "*2", new String[] { "C", "T", null, null });
    hap2.finalize(variants);

    NamedAllele hap3 = new NamedAllele("*3", "*3", new String[] { "C", null, null, "Y" });
    hap3.finalize(variants);

    // permutations should be consistent even if positions are out-of-order
    NamedAllele hap4 = new NamedAllele("*4", "*4", new String[] { "Y", null, null, "C" });
    Arrays.sort(variants, Collections.reverseOrder());
    hap4.finalize(variants);

    assertEquals("1:T;2:A;3:G;4:C;", ref.getPermutations().pattern());
    assertEquals("1:C;2:T;3:.?;4:.?;", hap2.getPermutations().pattern());
    assertEquals("1:C;2:.?;3:.?;4:[CT];", hap3.getPermutations().pattern());
    assertEquals("1:C;2:.?;3:.?;4:[CT];", hap4.getPermutations().pattern());
  }


  @Test
  public void testPermutationPattern() {

    String seq = "1:C;2:C;3:C;4:C;";
    Pattern p = Pattern.compile("1:C;2:.?;3:.?;4:[CT];");
    assertTrue(p.matcher(seq).matches());
  }
}

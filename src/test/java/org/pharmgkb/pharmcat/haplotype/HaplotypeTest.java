package org.pharmgkb.pharmcat.haplotype;

import java.util.Arrays;
import org.junit.Test;

import static org.junit.Assert.assertEquals;

/**
 * JUnit test for {@link Haplotype}.
 *
 * @author Mark Woon
 */
public class HaplotypeTest {


  @Test
  public void testCalculatePermutations() {

    Variant var1 = new Variant("chr1", "cyp2c19");
    var1.setPosition("1");

    Variant var2 = new Variant("chr1", "cyp2c19");
    var2.setPosition("2");

    Variant var3 = new Variant("chr1", "cyp2c19");
    var3.setPosition("3");

    Variant var4 = new Variant("chr1", "cyp2c19");
    var4.setPosition("4");

    Haplotype ref = new Haplotype(null, "*1");
    ref.addVariant(var1);
    ref.addVariant(var2);
    ref.addVariant(var3);
    ref.addVariant(var4);
    ref.addAlleles(Arrays.asList("T", "A", "G", "C"));

    Haplotype hap2 = new Haplotype(null, "*2");
    hap2.addVariant(var1);
    hap2.addVariant(var2);
    hap2.addAlleles(Arrays.asList("C", "T"));

    Haplotype hap3 = new Haplotype(null, "*3");
    hap3.addVariant(var1);
    hap3.addVariant(var4);
    hap3.addAlleles(Arrays.asList("C", "Y"));

    assertEquals("1:T;2:A;3:G;4:C;", ref.calculatePermutations(ref.getVariants()).pattern());
    assertEquals("1:C;2:T;3:.?;4:.?;", hap2.calculatePermutations(ref.getVariants()).pattern());
    assertEquals("1:C;2:.?;3:.?;4:C|T;", hap3.calculatePermutations(ref.getVariants()).pattern());
  }
}

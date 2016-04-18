package org.cpic.haplotype;

import java.util.Arrays;
import org.junit.Test;

import static org.junit.Assert.assertEquals;

/**
 * @author Mark Woon
 */
public class HaplotypeTest {


  @Test
  public void testCalculatePermutations() {

    Variant var1 = new Variant("chr1", "cyp2c19", "cdna");
    var1.setPOS("1");

    Variant var2 = new Variant("chr1", "cyp2c19", "cdna");
    var2.setPOS("2");

    Variant var3 = new Variant("chr1", "cyp2c19", "cdna");
    var3.setPOS("3");

    Variant var4 = new Variant("chr1", "cyp2c19", "cdna");
    var4.setPOS("4");

    Haplotype ref = new Haplotype();
    ref.setCommonName("*1");
    ref.addVariant(var1);
    ref.addVariant(var2);
    ref.addVariant(var3);
    ref.addVariant(var4);
    ref.addAlleles(Arrays.asList("T", "A", "G", "C"));

    Haplotype hap1 = new Haplotype();
    hap1.setCommonName("*2");
    hap1.addVariant(var1);
    hap1.addVariant(var2);
    hap1.addAlleles(Arrays.asList("C", "T"));

    Haplotype hap2 = new Haplotype();
    hap2.setCommonName("*2");
    hap2.addVariant(var1);
    hap2.addVariant(var4);
    hap2.addAlleles(Arrays.asList("C", "Y"));

    assertEquals("1:T;2:A;3:G;4:C;", ref.calculatePermutations(ref.getVariants()).pattern());
    assertEquals("1:C;2:T;3:.?;4:.?;", hap1.calculatePermutations(ref.getVariants()).pattern());
    assertEquals("1:C;2:.?;3:.?;4:C|T;", hap2.calculatePermutations(ref.getVariants()).pattern());
  }
}

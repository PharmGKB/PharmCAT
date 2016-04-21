package org.pharmgkb.pharmcat.haplotype;

import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.SortedSet;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import org.junit.Test;
import org.pharmgkb.pharmcat.haplotype.model.HaplotypeMatch;

import static org.junit.Assert.assertEquals;


/**
 * @author Mark Woon
 */
public class ComparisonUtilTest {


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

    SortedSet<HaplotypeMatch> matches = ComparisonUtil.comparePermutations(permutations, Lists.newArrayList(hap1, hap2, hap3));
    assertEquals(2, matches.size());
    Iterator<HaplotypeMatch> it = matches.iterator();
    assertEquals(hap1, it.next().getHaplotype());
    assertEquals(hap2, it.next().getHaplotype());
  }
}

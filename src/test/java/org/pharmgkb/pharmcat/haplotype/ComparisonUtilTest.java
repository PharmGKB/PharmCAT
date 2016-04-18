package org.pharmgkb.pharmcat.haplotype;

import java.util.List;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import org.junit.Test;


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

    Variant var1 = new Variant("chr1", "test", "cdna");
    var1.setPOS("1");

    Variant var2 = new Variant("chr1", "test", "cdna");
    var2.setPOS("2");

    Variant var3 = new Variant("chr1", "test", "cdna");
    var3.setPOS("3");

    Variant var4 = new Variant("chr1", "test", "cdna");
    var4.setPOS("4");

    List<Variant> variants = Lists.newArrayList(var1, var2, var3, var4);

    Haplotype hap1 = new Haplotype(
        variants,
        null,
        "*1",
        null,
        Lists.newArrayList("T", "A", "C", "C")
    );
    hap1.calculatePermutations(variants);

    Haplotype hap2 = new Haplotype(
        Lists.newArrayList(var2, var3),
        null,
        "*2",
        null,
        Lists.newArrayList("T", "C")
    );
    hap2.calculatePermutations(variants);

    Haplotype hap3 = new Haplotype(
        Lists.newArrayList(var4),
        null,
        "*3",
        null,
        Lists.newArrayList("GG")
    );
    hap3.calculatePermutations(variants);

    List<List<Haplotype>> pairs = CombinationUtil.generatePerfectPairs(new TreeSet<>(Sets.newHashSet(hap1, hap2, hap3)));

    SortedSet<HaplotypeMatch> matches = ComparisonUtil.comparePermutations(permutations, Lists.newArrayList(hap1, hap2, hap3));
    List<List<HaplotypeMatch>> pairMatches = ComparisonUtil.determinePairs(matches, pairs);
    ComparisonUtil.printMatchPairs(pairMatches);
  }
}

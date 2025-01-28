package org.pharmgkb.pharmcat.haplotype;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertSame;


/**
 * JUnit test for {@link SampleAllele}.
 *
 * @author Mark Woon
 */
class SampleAlleleTest {


  @Test
  void testCompare() {

    SampleAllele sa1 = new SampleAllele("chr1", 1, "A", "A", true, Lists.newArrayList("A", "G"), "0/0");
    SampleAllele sa2 = new SampleAllele("chr1", 2, "A", "A", true, Lists.newArrayList("A", "G"), "0/0");
    SampleAllele sa3 = new SampleAllele("chr2", 1, "A", "A", true, Lists.newArrayList("A", "G"), "0/0");
    SampleAllele sa4 = new SampleAllele("chrX", 1, "A", "A", true, Lists.newArrayList("A", "G"), "0/0");
    SampleAllele sa5 = new SampleAllele("chrY", 1, "A", "A", true, Lists.newArrayList("A", "G"), "0/0");
    SampleAllele sa6 = new SampleAllele("chr11", 1, "A", "A", true, Lists.newArrayList("A", "G"), "0/0");

    assertSame(Sets.newTreeSet(Sets.newHashSet(sa1, sa2)).first(), sa1);
    assertSame(Sets.newTreeSet(Sets.newHashSet(sa2, sa1)).first(), sa1);
    assertSame(Sets.newTreeSet(Sets.newHashSet(sa2, sa3)).first(), sa2);
    assertSame(Sets.newTreeSet(Sets.newHashSet(sa6, sa2)).first(), sa2);
    assertSame(Sets.newTreeSet(Sets.newHashSet(sa1, sa4)).first(), sa1);
    assertSame(Sets.newTreeSet(Sets.newHashSet(sa4, sa5)).first(), sa4);
    assertSame(Sets.newTreeSet(Sets.newHashSet(sa5, sa4)).first(), sa4);
  }
}

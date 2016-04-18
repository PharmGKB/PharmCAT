package org.pharmgkb.pharmcat;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import org.junit.Test;
import org.pharmgkb.pharmcat.haplotype.Haplotype;
import org.pharmgkb.pharmcat.haplotype.HaplotypeMatch;


/**
 * JUnit test for {@link TestUtil}.
 *
 * @author Mark Woon
 */
public class TestUtilTest {


  private HaplotypeMatch createHaplotypeMatch(String haplotypeName) {
    Haplotype hap = new Haplotype(null, null, haplotypeName, null, null);
    return new HaplotypeMatch(hap);
  }


  @Test
  @SuppressWarnings("unchecked")
  public void testAssertDiplotypePairs() {

    HaplotypeMatch hm1 = createHaplotypeMatch("*1");
    HaplotypeMatch hm2 = createHaplotypeMatch("*2");
    HaplotypeMatch hm3 = createHaplotypeMatch("*3");
    HaplotypeMatch hm4 = createHaplotypeMatch("*4");

    TestUtil.assertDiplotypePairs(
        Sets.newHashSet("*1/*1", "*2/*4"),
        Lists.newArrayList(
            Lists.newArrayList(hm1, hm1),
            Lists.newArrayList(hm2, hm4)
        )
    );

    TestUtil.assertDiplotypePairs(
        Sets.newHashSet("*1/*1", "*1/*3", "*2/*4", "*3/*4"),
        Lists.newArrayList(
            Lists.newArrayList(hm1, hm1),
            Lists.newArrayList(hm1, hm3),
            Lists.newArrayList(hm2, hm4),
            Lists.newArrayList(hm3, hm4)
        )
    );

    TestUtil.assertDiplotypePairs(
        Sets.newHashSet(),
        Lists.newArrayList(
        )
    );
  }

  @Test(expected = AssertionError.class)
  @SuppressWarnings("unchecked")
  public void testAssertDiplotypePairsMismatch() {

    HaplotypeMatch hm1 = createHaplotypeMatch("*1");
    HaplotypeMatch hm2 = createHaplotypeMatch("*2");

    TestUtil.assertDiplotypePairs(
        Sets.newHashSet("*1/*1", "*2/*4"),
        Lists.newArrayList(
            Lists.newArrayList(hm1, hm1),
            Lists.newArrayList(hm2, hm2)
        )
    );
  }
}

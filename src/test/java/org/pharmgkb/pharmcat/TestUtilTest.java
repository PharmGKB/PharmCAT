package org.pharmgkb.pharmcat;

import com.google.common.collect.Lists;
import org.junit.Test;
import org.pharmgkb.pharmcat.definition.model.NamedAllele;
import org.pharmgkb.pharmcat.haplotype.model.DiplotypeMatch;
import org.pharmgkb.pharmcat.haplotype.model.HaplotypeMatch;


/**
 * JUnit test for {@link TestUtil}.
 *
 * @author Mark Woon
 */
public class TestUtilTest {


  private HaplotypeMatch createHaplotypeMatch(String haplotypeName) {
    return new HaplotypeMatch(new NamedAllele(haplotypeName, haplotypeName, new String[1]));
  }


  @Test
  @SuppressWarnings("unchecked")
  public void testAssertDiplotypePairs() {

    HaplotypeMatch hm1 = createHaplotypeMatch("*1");
    HaplotypeMatch hm2 = createHaplotypeMatch("*2");
    HaplotypeMatch hm3 = createHaplotypeMatch("*3");
    HaplotypeMatch hm4 = createHaplotypeMatch("*4");

    TestUtil.assertDiplotypePairs(
        Lists.newArrayList("*1/*1", "*2/*4"),
        Lists.newArrayList(
            new DiplotypeMatch(hm1, hm1),
            new DiplotypeMatch(hm2, hm4)
        )
    );

    TestUtil.assertDiplotypePairs(
        Lists.newArrayList("*1/*1", "*1/*3", "*2/*4", "*3/*4"),
        Lists.newArrayList(
            new DiplotypeMatch(hm1, hm1),
            new DiplotypeMatch(hm1, hm3),
            new DiplotypeMatch(hm2, hm4),
            new DiplotypeMatch(hm3, hm4)
        )
    );

    TestUtil.assertDiplotypePairs(
        Lists.newArrayList(),
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
        Lists.newArrayList("*1/*1", "*2/*4"),
        Lists.newArrayList(
            new DiplotypeMatch(hm1, hm1),
            new DiplotypeMatch(hm2, hm2)
        )
    );
  }


  @Test
  @SuppressWarnings("unchecked")
  public void testAssertDiplotypePairsRightOrder() {

    HaplotypeMatch hm1 = createHaplotypeMatch("*1");
    HaplotypeMatch hm2 = createHaplotypeMatch("*2");

    TestUtil.assertDiplotypePairs(
        Lists.newArrayList("*2/*1", "*1/*2"),
        Lists.newArrayList(
            new DiplotypeMatch(hm2, hm1),
            new DiplotypeMatch(hm1, hm2)
        )
    );
  }

  @Test(expected = AssertionError.class)
  @SuppressWarnings("unchecked")
  public void testAssertDiplotypePairsWrongOrder() {

    HaplotypeMatch hm1 = createHaplotypeMatch("*1");
    HaplotypeMatch hm2 = createHaplotypeMatch("*2");

    TestUtil.assertDiplotypePairs(
        Lists.newArrayList("*1/*2", "*2/*1"),
        Lists.newArrayList(
            new DiplotypeMatch(hm2, hm1),
            new DiplotypeMatch(hm1, hm2)
        )
    );
  }
}

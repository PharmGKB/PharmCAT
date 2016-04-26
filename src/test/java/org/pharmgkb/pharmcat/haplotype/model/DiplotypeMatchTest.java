package org.pharmgkb.pharmcat.haplotype.model;

import java.util.SortedSet;
import java.util.TreeSet;
import com.google.common.collect.Lists;
import org.junit.Test;
import org.pharmgkb.pharmcat.haplotype.Haplotype;
import org.pharmgkb.pharmcat.haplotype.Variant;

import static org.junit.Assert.assertEquals;


/**
 * JUnit test for {@link DiplotypeMatch}.
 *
 * @author Mark Woon
 */
public class DiplotypeMatchTest {


  @Test
  public void testCompareTo() {

    Variant var1 = new Variant("chr1", "test");
    var1.setPosition("1");

    Variant var2 = new Variant("chr1", "test");
    var2.setPosition("2");

    Variant var3 = new Variant("chr1", "test");
    var3.setPosition("3");

    HaplotypeMatch hm1 = new HaplotypeMatch(new Haplotype(null, "*1", Lists.newArrayList(var1, var2, var3), null, null));
    HaplotypeMatch hm2 = new HaplotypeMatch(new Haplotype(null, "*4", Lists.newArrayList(var1, var2, var3), null, null));
    HaplotypeMatch hm3 = new HaplotypeMatch(new Haplotype(null, "*3", Lists.newArrayList(var1, var2), null, null));

    DiplotypeMatch dm1 = new DiplotypeMatch(hm1, hm1);
    DiplotypeMatch dm2 = new DiplotypeMatch(hm1, hm2);
    DiplotypeMatch dm3 = new DiplotypeMatch(hm2, hm2);
    DiplotypeMatch dm4 = new DiplotypeMatch(hm3, hm2);

    SortedSet<DiplotypeMatch> matches = new TreeSet<>(Lists.newArrayList(dm1, dm2));
    assertEquals(dm1, matches.first());

    matches = new TreeSet<>(Lists.newArrayList(dm3, dm2));
    assertEquals(dm2, matches.first());

    matches = new TreeSet<>(Lists.newArrayList(dm3, dm4));
    assertEquals(dm3, matches.first());
  }
}

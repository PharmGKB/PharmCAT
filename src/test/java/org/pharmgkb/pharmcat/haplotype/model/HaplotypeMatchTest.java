package org.pharmgkb.pharmcat.haplotype.model;

import java.util.SortedSet;
import java.util.TreeSet;
import org.junit.jupiter.api.Test;
import org.pharmgkb.pharmcat.definition.model.NamedAllele;

import static org.junit.jupiter.api.Assertions.assertEquals;


/**
 * JUnit test for {@link HaplotypeMatch}.
 *
 * @author Mark Woon
 */
class HaplotypeMatchTest {

  @Test
  void testCompare() {

    SortedSet<HaplotypeMatch> sortedSet = new TreeSet<>();
    HaplotypeMatch hm1 = new HaplotypeMatch(new NamedAllele("*4", "*4", new String[1]));
    sortedSet.add(hm1);
    HaplotypeMatch hm2 = new HaplotypeMatch(new NamedAllele("*1", "*1", new String[1]));
    sortedSet.add(hm2);
    assertEquals(hm2, sortedSet.first());
  }

}

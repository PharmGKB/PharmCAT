package org.pharmgkb.pharmcat.util;

import java.util.SortedSet;
import java.util.TreeSet;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;


/**
 * This is a JUnit test for {@link ChrNameComparator}.
 *
 * @author Mark Woon
 */
class ChrNameComparatorTest {

  @Test
  void testCompare() {
    ChrNameComparator comparator = new ChrNameComparator();

    SortedSet<String> chrs = new TreeSet<>(comparator);
    chrs.add("chr4");
    chrs.add("chrM");
    chrs.add("chr11");
    chrs.add("chrY");
    chrs.add("chr1");
    chrs.add("chrX");
    chrs.add("chr7");

    assertEquals("chr1, chr4, chr7, chr11, chrX, chrY, chrM", String.join(", ", chrs));
  }
}

package org.pharmgkb.pharmcat.util;

import java.util.SortedSet;
import java.util.TreeSet;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;


/**
 * This is a JUnit test for {@link HaplotypeNameComparator}.
 *
 * @author Mark Woon
 */
class HaplotypeNameComparatorTest {

  @Test
  void testBasic() {
    SortedSet<String> names = new TreeSet<>(new HaplotypeNameComparator());
    names.add("*3");
    names.add("*2");
    assertEquals("*2", names.first());
  }

  @Test
  void testComboComesLast() {
    SortedSet<String> names = new TreeSet<>(new HaplotypeNameComparator());
    names.add("[*1 + *2]");
    names.add("*2");
    assertEquals("*2", names.first());
  }

  @Test
  void testComboSortName() {
    SortedSet<String> names = new TreeSet<>(new HaplotypeNameComparator());
    names.add("[*3 + *4]");
    names.add("[*2 + *3]");
    assertEquals("[*2 + *3]", names.first());
  }

  @Test
  void testComboSortNumPartials() {
    SortedSet<String> names = new TreeSet<>(new HaplotypeNameComparator());
    names.add("[*2 + *3 + *4]");
    names.add("[*2 + *3]");
    assertEquals("[*2 + *3]", names.first());
  }
}

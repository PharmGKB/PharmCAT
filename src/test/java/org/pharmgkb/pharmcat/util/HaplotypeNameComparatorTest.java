package org.pharmgkb.pharmcat.util;

import java.util.SortedSet;
import java.util.TreeSet;
import org.junit.jupiter.api.Test;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.contains;
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


    names = new TreeSet<>(new HaplotypeNameComparator());
    names.add("[*2 + *4]");
    names.add("[*2 + *3]");
    assertEquals("[*2 + *3]", names.first());

    names = new TreeSet<>(new HaplotypeNameComparator());
    names.add("[c.85T>C (*9 A) + c.1218G>A + c.1627A>G (*5)]");
    names.add("[c.1218G>A + c.1627A>G (*5)]");
    System.out.println(names);
    assertThat(names, contains("[c.1218G>A + c.1627A>G (*5)]", "[c.85T>C (*9 A) + c.1218G>A + c.1627A>G (*5)]"));

    names = new TreeSet<>(new HaplotypeNameComparator());
    names.add("[c.1218G>A + c.1627A>G (*5)]");
    names.add("[*17 + c.1627A>G (*5)]");
    System.out.println(names);
    assertThat(names, contains("[*17 + c.1627A>G (*5)]", "[c.1218G>A + c.1627A>G (*5)]"));
  }

  @Test
  void testComboSortNumPartials() {
    SortedSet<String> names = new TreeSet<>(new HaplotypeNameComparator());
    names.add("[*2 + *3 + *4]");
    names.add("[*2 + *3]");
    assertEquals("[*2 + *3]", names.first());
  }

  @Test
  void testHgvsNames() {
    SortedSet<String> names = new TreeSet<>(new HaplotypeNameComparator());
    names.add("c.1627A>G (*5)");
    names.add("c.1371C>T");
    System.out.println(names);
    assertThat(names, contains("c.1371C>T", "c.1627A>G (*5)"));

    System.out.println();
    names.clear();
    names.add("c.1627A>G (*5)");
    names.add("c.85T>C (*9 A)");
    System.out.println(names);
    assertThat(names, contains("c.85T>C (*9 A)", "c.1627A>G (*5)"));
  }
}

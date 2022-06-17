package org.pharmgkb.pharmcat.reporter.model;

import java.util.Set;
import java.util.TreeSet;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;


class DrugLinkTest {

  @Test
  void testCompare() {
    DrugLink drugLink1 = new DrugLink("foo", "ID1");
    DrugLink drugLink2 = new DrugLink("bar", "ID2");
    DrugLink drugLink3 = new DrugLink("foo", "ID1");

    assertTrue(drugLink1.compareTo(drugLink2) > 0);
    assertEquals(0, drugLink1.compareTo(drugLink3));
    assertTrue(drugLink2.compareTo(drugLink3) < 0);

    Set<DrugLink> links = new TreeSet<>();
    links.add(drugLink1);
    links.add(drugLink2);
    links.add(drugLink3);
    assertEquals(2, links.size());
  }
}

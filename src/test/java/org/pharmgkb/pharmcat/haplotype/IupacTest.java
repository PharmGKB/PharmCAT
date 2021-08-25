package org.pharmgkb.pharmcat.haplotype;

import com.google.common.collect.Lists;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;


/**
 * This is a JUnit test for {@link Iupac}.
 *
 * @author Mark Woon
 */
class IupacTest {

  @Test
  void testGetBases() {
    assertEquals(Iupac.DEL.getBases(), Lists.newArrayList());
    assertEquals(Iupac.N.getBases(), Lists.newArrayList("A", "C", "G", "T"));
    assertEquals(Iupac.R.getBases(), Lists.newArrayList("A", "G"));
  }
}

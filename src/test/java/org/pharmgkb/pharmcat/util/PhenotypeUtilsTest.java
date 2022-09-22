package org.pharmgkb.pharmcat.util;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNull;


class PhenotypeUtilsTest {

  @Test
  void testNormalize() {
    assertEquals("Normal Metabolizer", PhenotypeUtils.normalize("Normal Metaboliser"));
    assertNull(PhenotypeUtils.normalize("  "));
    assertEquals("Normal Metabolizer", PhenotypeUtils.normalize("  Normal    Metabolizer"));
    assertEquals("Normal Metabolizer", PhenotypeUtils.normalize("NM"));
    assertEquals("Normal Metabolizer", PhenotypeUtils.normalize("EM"));
    assertEquals("Poor Metabolizer", PhenotypeUtils.normalize("PM"));
    assertEquals("Normal Metabolizer", PhenotypeUtils.normalize("Extensive Metaboliser"));
    assertEquals("Normal Metabolizer", PhenotypeUtils.normalize("Extensive Metabolizer"));
    assertEquals("Normal Metabolizer", PhenotypeUtils.normalize("Normal Metabolizer"));
    assertEquals("Normal Metabolizer", PhenotypeUtils.normalize("Normal MetaboliSer"));
    assertEquals("Normal Metabolizer", PhenotypeUtils.normalize("NM    "));
    assertEquals("Normal Metabolizer", PhenotypeUtils.normalize("normal metabolizer (NM)    "));
    assertEquals("Normal Metabolizer", PhenotypeUtils.normalize("normal   metabolizer (NM)\n"));
    assertEquals("foo bar", PhenotypeUtils.normalize("foo\t\t\t bar"));
    assertEquals("Likely Normal Metabolizer", PhenotypeUtils.normalize("Likely Normal Metaboliser"));
    assertEquals("Likely Normal Metabolizer", PhenotypeUtils.normalize("Likely Normal Metaboliser (NM)"));
  }
}

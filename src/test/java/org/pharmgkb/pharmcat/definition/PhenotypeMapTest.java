package org.pharmgkb.pharmcat.definition;

import org.junit.jupiter.api.Test;
import org.pharmgkb.pharmcat.definition.model.GenePhenotype;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;


/**
 * Test of the PhenotypeMap object to make sure it's loading properly
 *
 * @author Ryan Whaley
 */
class PhenotypeMapTest {

  @Test
  void test() throws Exception {
    PhenotypeMap phenotypeMap = new PhenotypeMap();

    assertNotNull(phenotypeMap);

    assertEquals(14, phenotypeMap.getGenes().size());

    assertEquals(
        "No function",
        phenotypeMap.lookup("CYP2C9").orElseThrow(Exception::new).getHaplotypes().get("*6"));

    GenePhenotype genePhenotype = phenotypeMap.lookup("DPYD").orElseThrow(Exception::new);
    assertNotNull(genePhenotype);
    assertEquals("Normal function", genePhenotype.lookupHaplotype("Reference"));
  }

  @Test
  void testLookupPhenotype() throws Exception {
    PhenotypeMap phenotypeMap = new PhenotypeMap();

    GenePhenotype genePhenotype = phenotypeMap.lookup("CYP2C9")
        .orElseThrow(() -> new RuntimeException("No CYP2C9 phenotype map found"));
    assertNotNull(genePhenotype.getDiplotypes());
    assertEquals("2", genePhenotype.getLookupKeyForDiplotype("*1/*1"));
  }
}

package org.pharmgkb.pharmcat.definition;

import org.junit.Test;
import org.pharmgkb.pharmcat.definition.model.GenePhenotype;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;


/**
 * Test of the PhenotypeMap object to make sure it's loading properly
 *
 * @author Ryan Whaley
 */
public class PhenotypeMapTest {

  @Test
  public void test() throws Exception {
    PhenotypeMap phenotypeMap = new PhenotypeMap();

    assertNotNull(phenotypeMap);

    assertEquals(12, phenotypeMap.getGenes().size());

    assertEquals(
        "No function",
        phenotypeMap.lookup("CYP2C9").orElseThrow(Exception::new).getHaplotypes().get("*6"));

    GenePhenotype genePhenotype = phenotypeMap.lookup("DPYD").orElseThrow(Exception::new);
    assertNotNull(genePhenotype);
    assertEquals("Normal function", genePhenotype.lookupHaplotype("Any normal function variant or no variant detected"));
  }
}

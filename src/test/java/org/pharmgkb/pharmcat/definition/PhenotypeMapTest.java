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
    PhenotypeMap phenotypeMap = new PhenotypeMap(null);

    assertNotNull(phenotypeMap);

    assertEquals(8, phenotypeMap.getGenes().size());

    assertEquals(
        "No Function",
        phenotypeMap.lookup("CYP2C9").orElseThrow(Exception::new).getHaplotypes().get("*6"));

    GenePhenotype genePhenotype = phenotypeMap.lookup("DPYD").orElseThrow(Exception::new);
    assertNotNull(genePhenotype);
    assertEquals("Normal Function", genePhenotype.lookupHaplotype("Any normal function variant or no variant detected"));
  }
}

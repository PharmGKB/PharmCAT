package org.pharmgkb.pharmcat.definition;

import org.junit.Test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;


/**
 * Test of the PhenotypeMap object to make sure it's loading properly
 *
 * @author whaleyr
 */
public class PhenotypeMapTest {

  @Test
  public void test() throws Exception {
    PhenotypeMap phenotypeMap = new PhenotypeMap();

    assertNotNull(phenotypeMap);

    assertEquals(8, phenotypeMap.getGenes().size());

    assertEquals("No Function", phenotypeMap.lookup("CYP2C9").getHaplotypes().get("*6"));
  }
}

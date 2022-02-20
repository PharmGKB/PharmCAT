package org.pharmgkb.pharmcat.definition;

import org.junit.jupiter.api.Test;
import org.pharmgkb.pharmcat.definition.model.GenePhenotype;
import org.pharmgkb.pharmcat.reporter.model.result.Diplotype;
import org.pharmgkb.pharmcat.reporter.model.result.Haplotype;

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

    assertEquals(16, phenotypeMap.getGenes().size());

    assertEquals(
        "No function",
        phenotypeMap.lookup("CYP2C9").orElseThrow(Exception::new).getHaplotypes().get("*6"));

    GenePhenotype genePhenotype = phenotypeMap.lookup("DPYD").orElseThrow(Exception::new);
    assertNotNull(genePhenotype);
    assertEquals("Normal function", genePhenotype.lookupHaplotype("Reference"));
  }

  @Test
  void testLookupPhenotype() {
    PhenotypeMap phenotypeMap = new PhenotypeMap();
    Diplotype diplotype = new Diplotype(
        "CYP2C19",
        new Haplotype("CYP2C19", "*1"),
        new Haplotype("CYP2C19", "*1")
    );

    GenePhenotype genePhenotype = phenotypeMap.lookup("CYP2C9")
        .orElseThrow(() -> new RuntimeException("No CYP2C9 phenotype map found"));
    assertNotNull(genePhenotype.getDiplotypes());
    assertEquals("2", genePhenotype.getLookupKeyForDiplotype(diplotype));
  }

  @Test
  void testLookupDpyd() {
    PhenotypeMap phenotypeMap = new PhenotypeMap();

    GenePhenotype genePhenotype = phenotypeMap.lookup("DPYD")
        .orElseThrow(() -> new RuntimeException("No DPYD phenotype map found"));
    assertNotNull(genePhenotype.getDiplotypes());

    Diplotype diplotype1 = new Diplotype(
        "DPYD",
        new Haplotype("DPYD", "Reference"),
        new Haplotype("DPYD", "c.2846A>T")
    );
    assertEquals("1.5", genePhenotype.getLookupKeyForDiplotype(diplotype1));

    Diplotype diplotype2 = new Diplotype(
        "DPYD",
        new Haplotype("DPYD", "c.2846A>T"),
        new Haplotype("DPYD", "Reference")
    );
    assertEquals("1.5", genePhenotype.getLookupKeyForDiplotype(diplotype2));

    Diplotype diplotype3 = new Diplotype(
        "DPYD",
        new Haplotype("DPYD", "foo"),
        new Haplotype("DPYD", "bar")
    );
    assertEquals("N/A", genePhenotype.getLookupKeyForDiplotype(diplotype3));
  }
}

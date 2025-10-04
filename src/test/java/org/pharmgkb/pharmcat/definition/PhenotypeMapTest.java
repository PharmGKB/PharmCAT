package org.pharmgkb.pharmcat.definition;

import java.util.Collection;
import java.util.List;
import java.util.Objects;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import org.pharmgkb.pharmcat.Constants;
import org.pharmgkb.pharmcat.Env;
import org.pharmgkb.pharmcat.phenotype.PhenotypeMap;
import org.pharmgkb.pharmcat.phenotype.model.DiplotypeRecord;
import org.pharmgkb.pharmcat.phenotype.model.GenePhenotype;
import org.pharmgkb.pharmcat.reporter.TextConstants;
import org.pharmgkb.pharmcat.reporter.model.DataSource;
import org.pharmgkb.pharmcat.reporter.model.result.Diplotype;
import org.pharmgkb.pharmcat.reporter.model.result.Haplotype;

import static org.junit.jupiter.api.Assertions.*;


/**
 * Test of the PhenotypeMap object to make sure it's loading properly
 *
 * @author Ryan Whaley
 */
class PhenotypeMapTest {
  private static Env s_env;

  @BeforeAll
  static void prepare() throws Exception {
    s_env = new Env();
  }


  @Test
  void test() throws Exception {
    PhenotypeMap phenotypeMap = new PhenotypeMap();

    assertNotNull(phenotypeMap);
    List<String> genes = List.of(
        "ABCG2",
        "CACNA1S",
        "CFTR",
        "CYP2B6",
        "CYP2C19",
        "CYP2C9",
        "CYP2D6",
        "CYP3A5",
        "DPYD",
        "G6PD",
        "IFNL3",
        "MT-RNR1",
        "NAT2",
        "NUDT15",
        "RYR1",
        "SLCO1B1",
        "TPMT",
        "UGT1A1"
    );
    List<String> cpicGenes = phenotypeMap.getCpicGenes().stream()
        .map(GenePhenotype::getGene)
        .sorted()
        .toList();
    assertEquals(genes, cpicGenes);

    // HLA's are not part of the Phenotype map, they use allele status instead
    assertTrue(phenotypeMap.getCpicGenes().stream().noneMatch((gene) -> gene.getGene().startsWith("HLA")));

    assertEquals(
        "0.0",
        phenotypeMap.lookupCpic("CYP2C9").orElseThrow(Exception::new).getHaplotypes().get("*6"));

    GenePhenotype genePhenotype = phenotypeMap.lookupCpic("DPYD").orElseThrow(Exception::new);
    assertNotNull(genePhenotype);
    assertEquals("Normal function", genePhenotype.getHaplotypeFunction("Reference"));
  }

  /**
   * Ensure all diplotype calculation gets assigned lookup keys where appropriate
   */
  @Test
  void testLookupTranslation() {
    DataSource[] sources = new DataSource[] {DataSource.CPIC, DataSource.DPWG};
    PhenotypeMap phenotypeMap = new PhenotypeMap();

    for (DataSource source : sources) {
      Collection<GenePhenotype> allGenePhenotypes = source == DataSource.CPIC ? phenotypeMap.getCpicGenes() : phenotypeMap.getDpwgGenes();
      for (GenePhenotype genePhenotype : allGenePhenotypes) {
        boolean isScoreGene = Constants.isActivityScoreGene(genePhenotype.getGene(), source);
        String unspecifiedString = isScoreGene ? "n/a" : "Indeterminate";
        for (DiplotypeRecord diplotype : genePhenotype.getDiplotypes()) {
          String description = source + " " + genePhenotype + " " + diplotype;
          if (diplotype.getPhenotype().equals("Indeterminate")) {
            assertEquals("Indeterminate", diplotype.getGeneResult(), "No gene result not set to Indeterminate for " + description);
            assertEquals(unspecifiedString, diplotype.getLookupKey(), "Lookup key should be n/a for " + description);
            assertTrue(!isScoreGene || TextConstants.isUnspecified(diplotype.getActivityScore()), "Activity score should be n/a for " + description);
          } else {
            assertFalse(TextConstants.isUnspecified(diplotype.getPhenotype()), "No phenotype for " + description);
            assertFalse(TextConstants.isUnspecified(diplotype.getGeneResult()), "No gene result for " + description);
            assertFalse(TextConstants.isUnspecified(diplotype.getLookupKey()), "No lookup key for " + description);
            assertFalse(isScoreGene && TextConstants.isUnspecified(diplotype.getActivityScore()), "No activity score for " + description);
          }
        }
      }
    }
  }

  @Test
  void testLookupPhenotype() {
    PhenotypeMap phenotypeMap = new PhenotypeMap();
    Diplotype diplotype = new Diplotype("CYP2C19", "*1", "*1", s_env, DataSource.CPIC, 0);

    GenePhenotype genePhenotype = phenotypeMap.lookupCpic("CYP2C19")
        .orElseThrow(() -> new RuntimeException("No CYP2C19 phenotype map found"));
    assertNotNull(genePhenotype.getDiplotypes());

    assertEquals(1, diplotype.getLookupKeys().size());
    assertEquals("Normal Metabolizer", diplotype.getLookupKeys().get(0));
  }

  @Test
  void testLookupActivity() {
    PhenotypeMap phenotypeMap = new PhenotypeMap();
    Diplotype diplotype = new Diplotype("CYP2D6", "*1", "*3", s_env, DataSource.CPIC, 0);

    GenePhenotype genePhenotype = phenotypeMap.lookupCpic("CYP2D6")
        .orElseThrow(() -> new RuntimeException("No CYP2D6 phenotype map found"));
    assertNotNull(genePhenotype.getDiplotypes());

    assertTrue(genePhenotype.getDiplotypes().stream()
        .anyMatch((d) -> Objects.nonNull(d.getActivityScore())));

    assertEquals("1.0", diplotype.getActivityScore());
  }

  @Test
  void testLookupDpyd() {
    PhenotypeMap phenotypeMap = new PhenotypeMap();

    GenePhenotype genePhenotype = phenotypeMap.lookupCpic("DPYD")
        .orElseThrow(() -> new RuntimeException("No DPYD phenotype map found"));
    assertNotNull(genePhenotype.getDiplotypes());

    Diplotype diplotype1 = new Diplotype(
        "DPYD",
        new Haplotype("DPYD", "Reference"),
        new Haplotype("DPYD", "c.2846A>T"),
        s_env, DataSource.CPIC
    );
    assertEquals(1, diplotype1.getLookupKeys().size());
    assertEquals("1.5", diplotype1.getLookupKeys().get(0));

    Diplotype diplotype2 = new Diplotype(
        "DPYD",
        new Haplotype("DPYD", "c.2846A>T"),
        new Haplotype("DPYD", "Reference"),
        s_env, DataSource.CPIC
    );
    assertEquals(1, diplotype2.getLookupKeys().size());
    assertEquals("1.5", diplotype2.getLookupKeys().get(0));

    Diplotype diplotype3 = new Diplotype(
        "DPYD",
        new Haplotype("DPYD", "foo"),
        new Haplotype("DPYD", "bar"),
        s_env, DataSource.CPIC
    );
    assertEquals(1, diplotype3.getLookupKeys().size());
    assertEquals(TextConstants.NA, diplotype3.getLookupKeys().get(0));
  }
}

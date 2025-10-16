package org.pharmgkb.pharmcat.reporter.model.result;

import java.util.List;
import com.google.common.collect.ImmutableList;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import org.pharmgkb.pharmcat.Env;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;


public class GenotypeTest {
  private static final String sf_gene1 = "CYP2D6";
  private static final Haplotype sf_gene1h1 = new Haplotype(sf_gene1, "*1");
  private static final Haplotype sf_gene1h2 = new Haplotype(sf_gene1, "*3");

  private static final String sf_gene2 = "GENEX";
  private static final Haplotype sf_gene2h1 = new Haplotype(sf_gene2, "*1");
  private static final Haplotype sf_gene2h2 = new Haplotype(sf_gene2, "*2");

  private static Diplotype s_diplotype1;
  private static Diplotype s_diplotype2;
  private static Diplotype s_diplotype3;


  @BeforeAll
  static void prepare() throws Exception {
    Env env = new Env();
    s_diplotype1 = new Diplotype(sf_gene1, sf_gene1h1, sf_gene1h2, env);
    s_diplotype2 = new Diplotype(sf_gene1, sf_gene1h1, sf_gene1h1, env);
    s_diplotype3 = new Diplotype(sf_gene2, sf_gene2h1, sf_gene2h2, env);
  }




  @Test
  void testMakeGenotypes() {
    GeneReport geneReport1 = new GeneReport(sf_gene1, "test");
    geneReport1.addReporterDiplotype(s_diplotype1);
    geneReport1.addReporterDiplotype(s_diplotype2);
    GeneReport geneReport2 = new GeneReport(sf_gene2, "test");
    geneReport2.addReporterDiplotype(s_diplotype3);
    List<GeneReport> geneReports = ImmutableList.of(geneReport1, geneReport2);
    List<Genotype> possibleGenotypes = Genotype.makeGenotypes(geneReports);

    assertNotNull(possibleGenotypes);
    assertEquals(2, possibleGenotypes.size());

    for (Genotype possibleGenotype : possibleGenotypes) {
      assertEquals(1, possibleGenotype.getDiplotypes().stream().filter(d -> d.getGene().equals(sf_gene1)).count());
      assertEquals(1, possibleGenotype.getDiplotypes().stream().filter(d -> d.getGene().equals(sf_gene2)).count());
    }
  }

  @Test
  void testMakeGenotypesOneDiplotype() {
    GeneReport geneReport1 = new GeneReport(sf_gene1, "test");
    geneReport1.addReporterDiplotype(s_diplotype1);
    List<GeneReport> geneReports = ImmutableList.of(geneReport1);
    List<Genotype> possibleGenotypes = Genotype.makeGenotypes(geneReports);

    assertNotNull(possibleGenotypes);
    assertEquals(1, possibleGenotypes.size());

    for (Genotype possibleGenotype : possibleGenotypes) {
      assertEquals(1, possibleGenotype.getDiplotypes().stream().filter(d -> d.getGene().equals(sf_gene1)).count());
      assertEquals(0, possibleGenotype.getDiplotypes().stream().filter(d -> d.getGene().equals(sf_gene2)).count());
    }
  }

  @Test
  void testMakeGenotypesOneGeneTwoDiplotypes() {
    GeneReport geneReport1 = new GeneReport(sf_gene1, "test");
    geneReport1.addReporterDiplotype(s_diplotype1);
    geneReport1.addReporterDiplotype(s_diplotype2);
    List<GeneReport> geneReports = ImmutableList.of(geneReport1);
    List<Genotype> possibleGenotypes = Genotype.makeGenotypes(geneReports);

    assertNotNull(possibleGenotypes);
    assertEquals(2, possibleGenotypes.size());

    for (Genotype possibleGenotype : possibleGenotypes) {
      assertEquals(1, possibleGenotype.getDiplotypes().stream().filter(d -> d.getGene().equals(sf_gene1)).count());
      assertEquals(0, possibleGenotype.getDiplotypes().stream().filter(d -> d.getGene().equals(sf_gene2)).count());
    }
  }

  @Test
  void testMakeGenotypesTwoGenes() {
    GeneReport geneReport1 = new GeneReport(sf_gene1, "test");
    geneReport1.addReporterDiplotype(s_diplotype1);
    GeneReport geneReport2 = new GeneReport(sf_gene2, "test");
    geneReport2.addReporterDiplotype(s_diplotype3);
    List<GeneReport> geneReports = ImmutableList.of(geneReport1, geneReport2);
    List<Genotype> possibleGenotypes = Genotype.makeGenotypes(geneReports);

    assertNotNull(possibleGenotypes);
    assertEquals(1, possibleGenotypes.size());

    for (Genotype possibleGenotype : possibleGenotypes) {
      assertEquals(1, possibleGenotype.getDiplotypes().stream().filter(d -> d.getGene().equals(sf_gene1)).count());
      assertEquals(1, possibleGenotype.getDiplotypes().stream().filter(d -> d.getGene().equals(sf_gene2)).count());
    }
  }
}

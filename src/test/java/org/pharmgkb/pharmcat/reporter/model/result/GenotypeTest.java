package org.pharmgkb.pharmcat.reporter.model.result;

import java.util.List;
import com.google.common.collect.ImmutableList;
import org.junit.jupiter.api.Test;
import org.pharmgkb.pharmcat.reporter.model.DataSource;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;


public class GenotypeTest {
  private static final String gene1 = "CYP2D6";
  private static final Haplotype gene1h1 = new Haplotype(gene1, "*1");
  private static final Haplotype gene1h2 = new Haplotype(gene1, "*3");
  private static final Diplotype diplotype1 = new Diplotype(gene1, gene1h1, gene1h2);
  private static final Diplotype diplotype2 = new Diplotype(gene1, gene1h1, gene1h1);

  private static final String gene2 = "GENEX";
  private static final Haplotype gene2h1 = new Haplotype(gene2, "*1");
  private static final Haplotype gene2h2 = new Haplotype(gene2, "*2");
  private static final Diplotype diplotype3 = new Diplotype(gene2, gene2h1, gene2h2);


  @Test
  void testMakeGenotypes() {
    GeneReport geneReport1 = new GeneReport(gene1, DataSource.CPIC, "test");
    geneReport1.addReporterDiplotype(diplotype1);
    geneReport1.addReporterDiplotype(diplotype2);
    GeneReport geneReport2 = new GeneReport(gene2, DataSource.CPIC, "test");
    geneReport2.addReporterDiplotype(diplotype3);
    List<GeneReport> geneReports = ImmutableList.of(geneReport1, geneReport2);
    List<Genotype> possibleGenotypes = Genotype.makeGenotypes(geneReports);

    assertNotNull(possibleGenotypes);
    assertEquals(2, possibleGenotypes.size());

    for (Genotype possibleGenotype : possibleGenotypes) {
      assertEquals(1, possibleGenotype.getDiplotypes().stream().filter(d -> d.getGene().equals(gene1)).count());
      assertEquals(1, possibleGenotype.getDiplotypes().stream().filter(d -> d.getGene().equals(gene2)).count());
    }
  }

  @Test
  void testMakeGenotypesOneDiplotype() {
    GeneReport geneReport1 = new GeneReport(gene1, DataSource.CPIC, "test");
    geneReport1.addReporterDiplotype(diplotype1);
    List<GeneReport> geneReports = ImmutableList.of(geneReport1);
    List<Genotype> possibleGenotypes = Genotype.makeGenotypes(geneReports);

    assertNotNull(possibleGenotypes);
    assertEquals(1, possibleGenotypes.size());

    for (Genotype possibleGenotype : possibleGenotypes) {
      assertEquals(1, possibleGenotype.getDiplotypes().stream().filter(d -> d.getGene().equals(gene1)).count());
      assertEquals(0, possibleGenotype.getDiplotypes().stream().filter(d -> d.getGene().equals(gene2)).count());
    }
  }

  @Test
  void testMakeGenotypesOneGeneTwoDiplotypes() {
    GeneReport geneReport1 = new GeneReport(gene1, DataSource.CPIC, "test");
    geneReport1.addReporterDiplotype(diplotype1);
    geneReport1.addReporterDiplotype(diplotype2);
    List<GeneReport> geneReports = ImmutableList.of(geneReport1);
    List<Genotype> possibleGenotypes = Genotype.makeGenotypes(geneReports);

    assertNotNull(possibleGenotypes);
    assertEquals(2, possibleGenotypes.size());

    for (Genotype possibleGenotype : possibleGenotypes) {
      assertEquals(1, possibleGenotype.getDiplotypes().stream().filter(d -> d.getGene().equals(gene1)).count());
      assertEquals(0, possibleGenotype.getDiplotypes().stream().filter(d -> d.getGene().equals(gene2)).count());
    }
  }

  @Test
  void testMakeGenotypesTwoGenes() {
    GeneReport geneReport1 = new GeneReport(gene1, DataSource.CPIC, "test");
    geneReport1.addReporterDiplotype(diplotype1);
    GeneReport geneReport2 = new GeneReport(gene2, DataSource.CPIC, "test");
    geneReport2.addReporterDiplotype(diplotype3);
    List<GeneReport> geneReports = ImmutableList.of(geneReport1, geneReport2);
    List<Genotype> possibleGenotypes = Genotype.makeGenotypes(geneReports);

    assertNotNull(possibleGenotypes);
    assertEquals(1, possibleGenotypes.size());

    for (Genotype possibleGenotype : possibleGenotypes) {
      assertEquals(1, possibleGenotype.getDiplotypes().stream().filter(d -> d.getGene().equals(gene1)).count());
      assertEquals(1, possibleGenotype.getDiplotypes().stream().filter(d -> d.getGene().equals(gene2)).count());
    }
  }
}

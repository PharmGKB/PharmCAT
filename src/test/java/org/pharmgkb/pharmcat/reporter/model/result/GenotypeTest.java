package org.pharmgkb.pharmcat.reporter.model.result;

import java.util.List;
import com.google.common.collect.ImmutableList;
import org.junit.jupiter.api.Test;

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
    List<Diplotype> possibleDiplotypes = ImmutableList.of(diplotype1, diplotype2, diplotype3);
    List<Genotype> possibleGenotypes = Genotype.makeGenotypes(possibleDiplotypes);

    assertNotNull(possibleGenotypes);
    assertEquals(2, possibleGenotypes.size());

    for (Genotype possibleGenotype : possibleGenotypes) {
      assertEquals(1, possibleGenotype.getDiplotypes().stream().filter(d -> d.getGene().equals(gene1)).count());
      assertEquals(1, possibleGenotype.getDiplotypes().stream().filter(d -> d.getGene().equals(gene2)).count());
    }
  }

  @Test
  void testMakeGenotypesOneDiplotype() {
    List<Diplotype> possibleDiplotypes = ImmutableList.of(diplotype1);
    List<Genotype> possibleGenotypes = Genotype.makeGenotypes(possibleDiplotypes);

    assertNotNull(possibleGenotypes);
    assertEquals(1, possibleGenotypes.size());

    for (Genotype possibleGenotype : possibleGenotypes) {
      assertEquals(1, possibleGenotype.getDiplotypes().stream().filter(d -> d.getGene().equals(gene1)).count());
      assertEquals(0, possibleGenotype.getDiplotypes().stream().filter(d -> d.getGene().equals(gene2)).count());
    }
  }

  @Test
  void testMakeGenotypesOneGeneTwoDiplotypes() {
    List<Diplotype> possibleDiplotypes = ImmutableList.of(diplotype1, diplotype2);
    List<Genotype> possibleGenotypes = Genotype.makeGenotypes(possibleDiplotypes);

    assertNotNull(possibleGenotypes);
    assertEquals(2, possibleGenotypes.size());

    for (Genotype possibleGenotype : possibleGenotypes) {
      assertEquals(1, possibleGenotype.getDiplotypes().stream().filter(d -> d.getGene().equals(gene1)).count());
      assertEquals(0, possibleGenotype.getDiplotypes().stream().filter(d -> d.getGene().equals(gene2)).count());
    }
  }

  @Test
  void testMakeGenotypesTwoGenes() {
    List<Diplotype> possibleDiplotypes = ImmutableList.of(diplotype1, diplotype3);
    List<Genotype> possibleGenotypes = Genotype.makeGenotypes(possibleDiplotypes);

    assertNotNull(possibleGenotypes);
    assertEquals(1, possibleGenotypes.size());

    for (Genotype possibleGenotype : possibleGenotypes) {
      assertEquals(1, possibleGenotype.getDiplotypes().stream().filter(d -> d.getGene().equals(gene1)).count());
      assertEquals(1, possibleGenotype.getDiplotypes().stream().filter(d -> d.getGene().equals(gene2)).count());
    }
  }
}

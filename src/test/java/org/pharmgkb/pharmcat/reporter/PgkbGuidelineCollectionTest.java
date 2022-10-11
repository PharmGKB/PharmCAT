package org.pharmgkb.pharmcat.reporter;

import java.io.IOException;
import java.util.Set;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import org.pharmgkb.pharmcat.Env;
import org.pharmgkb.pharmcat.reporter.model.DataSource;
import org.pharmgkb.pharmcat.reporter.model.result.Diplotype;
import org.pharmgkb.pharmcat.reporter.model.result.Haplotype;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;


class PgkbGuidelineCollectionTest {
  private static Env s_env;

  @BeforeAll
  static void prepare() throws Exception {
    s_env = new Env();
  }

  @Test
  void testLoad() throws IOException {
    PgkbGuidelineCollection pgkbGuidelineCollection = new PgkbGuidelineCollection();
    System.out.println("# DPWG guideilnes: " + pgkbGuidelineCollection.getGuidelinePackages().size());
    assertTrue(pgkbGuidelineCollection.getGuidelinePackages().size() > 50);
  }

  @Test
  void testGetPhenotypesForDiplotype() throws IOException {
    String gene = "CYP2C19";
    Diplotype diplotype = new Diplotype(gene, new Haplotype(gene, "*1"), new Haplotype(gene, "*1"), s_env,
        DataSource.CPIC);

    PgkbGuidelineCollection pgkbGuidelineCollection = new PgkbGuidelineCollection();
    Set<String> phenotypes = pgkbGuidelineCollection.getPhenotypesForDiplotype(diplotype);

    assertEquals(1, phenotypes.size());

    assertEquals("Normal Metabolizer", String.join("::", phenotypes));
  }
}

package org.pharmgkb.pharmcat.reporter;

import java.io.IOException;
import java.util.Set;

import com.google.common.html.HtmlEscapers;
import org.junit.jupiter.api.Test;
import org.pharmgkb.pharmcat.reporter.model.result.Diplotype;
import org.pharmgkb.pharmcat.reporter.model.result.Haplotype;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;


class PgkbGuidelineCollectionTest {

  @Test
  void testLoad() throws IOException {
    PgkbGuidelineCollection pgkbGuidelineCollection = new PgkbGuidelineCollection();
    assertEquals(62, pgkbGuidelineCollection.getGuidelinePackages().size());
  }

  @Test
  void testGetPhenotypesForDiplotype() throws IOException {
    String gene = "CYP2C19";
    Diplotype diplotype = new Diplotype(gene, new Haplotype(gene, "*1"), new Haplotype(gene, "*1"));

    PgkbGuidelineCollection pgkbGuidelineCollection = new PgkbGuidelineCollection();
    Set<String> phenotypes = pgkbGuidelineCollection.getPhenotypesForDiplotype(diplotype);

    assertEquals(1, phenotypes.size());

    assertEquals("Normal Metabolizer", String.join("::", phenotypes));
  }
}

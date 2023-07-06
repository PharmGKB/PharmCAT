package org.pharmgkb.pharmcat.reporter;

import java.io.IOException;
import java.util.List;
import java.util.Set;
import java.util.SortedSet;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import org.pharmgkb.pharmcat.Env;
import org.pharmgkb.pharmcat.reporter.model.DataSource;
import org.pharmgkb.pharmcat.reporter.model.pgkb.GuidelinePackage;
import org.pharmgkb.pharmcat.reporter.model.result.Diplotype;
import org.pharmgkb.pharmcat.reporter.model.result.Haplotype;

import static org.junit.jupiter.api.Assertions.*;


class PgkbGuidelineCollectionTest {
  private static Env s_env;

  @BeforeAll
  static void prepare() throws Exception {
    s_env = new Env();
  }

  @Test
  void testLoad() throws IOException {
    PgkbGuidelineCollection pgkbGuidelineCollection = new PgkbGuidelineCollection();
    int guidelineCount = pgkbGuidelineCollection.getGuidelinePackages().size();
    System.out.println("# PharmGKB guideilnes: " + guidelineCount);
    assertTrue(guidelineCount > 50);

    Set<GuidelinePackage> cpicGuidelines = pgkbGuidelineCollection.getGuidelinesFromSource(DataSource.CPIC);
    Set<GuidelinePackage> dpwgGuidelines = pgkbGuidelineCollection.getGuidelinesFromSource(DataSource.DPWG);

    assertEquals(cpicGuidelines.size() + dpwgGuidelines.size(), guidelineCount);

    SortedSet<String> dpwgGenes = pgkbGuidelineCollection.getGenesUsedInSource(DataSource.DPWG);
    assertFalse(dpwgGenes.contains("CACNA1S"));
  }
}

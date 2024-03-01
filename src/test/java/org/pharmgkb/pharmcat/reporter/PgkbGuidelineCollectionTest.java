package org.pharmgkb.pharmcat.reporter;

import java.io.IOException;
import java.util.Set;
import java.util.SortedSet;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import org.pharmgkb.pharmcat.Env;
import org.pharmgkb.pharmcat.reporter.model.DataSource;
import org.pharmgkb.pharmcat.reporter.model.PrescribingGuidanceSource;
import org.pharmgkb.pharmcat.reporter.model.pgkb.GuidelinePackage;

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
    System.out.println("# PharmGKB prescribing guidance documents: " + guidelineCount);
    assertTrue(guidelineCount > 50);

    Set<GuidelinePackage> cpicGuidelines = pgkbGuidelineCollection.getGuidelinesFromSource(PrescribingGuidanceSource.CPIC_GUIDELINE);
    Set<GuidelinePackage> dpwgGuidelines = pgkbGuidelineCollection.getGuidelinesFromSource(PrescribingGuidanceSource.DPWG_GUIDELINE);
    Set<GuidelinePackage> fdaLabels = pgkbGuidelineCollection.getGuidelinesFromSource(PrescribingGuidanceSource.FDA_LABEL);
    Set<GuidelinePackage> fdaAssocs = pgkbGuidelineCollection.getGuidelinesFromSource(PrescribingGuidanceSource.FDA_ASSOC);

    assertEquals(cpicGuidelines.size() + dpwgGuidelines.size() + fdaLabels.size() + fdaAssocs.size(), guidelineCount);

    SortedSet<String> dpwgGenes = pgkbGuidelineCollection.getGenesUsedInSource(DataSource.DPWG);
    assertFalse(dpwgGenes.contains("CACNA1S"));
  }
}

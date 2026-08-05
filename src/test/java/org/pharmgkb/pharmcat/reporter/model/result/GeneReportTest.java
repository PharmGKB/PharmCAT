package org.pharmgkb.pharmcat.reporter.model.result;

import java.util.List;
import java.util.TreeSet;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import org.pharmgkb.pharmcat.Env;
import org.pharmgkb.pharmcat.phenotype.model.OutsideCall;
import org.pharmgkb.pharmcat.reporter.format.html.ReportHelpers;

import static org.junit.jupiter.api.Assertions.*;

class GeneReportTest {
  private static Env s_env;

  @BeforeAll
  static void prepare() throws Exception {
    s_env = new Env();
  }


  @Test
  void testPlaceholderGene() {
    String gene = "GENEX";
    GeneReport geneReport = new GeneReport(gene, "test");
    assertEquals(gene, geneReport.getGene());
    assertFalse(geneReport.isCalled());
    assertEquals(0, geneReport.getRecommendationDiplotypes().size());
  }

  @Test
  void testOutsideGene() {
    String gene = "UGT1A1";
    String outsideCallData = "UGT1A1\t*1/*6";
    String displayDiplotype = "*1/*6";
    String diplotypeString = "UGT1A1:*1/*6";

    OutsideCall outsideCall = new OutsideCall(s_env, outsideCallData, 0);
    GeneReport geneReport = new GeneReport(outsideCall, s_env);

    assertEquals(gene, geneReport.getGene());
    assertTrue(geneReport.isReportable());
    assertEquals(1, geneReport.getRecommendationDiplotypes().size());
    List<String> geneCalls = ReportHelpers.amdGeneCalls(geneReport);
    assertEquals(1, geneCalls.size());
    assertTrue(geneCalls.contains(displayDiplotype));

    Diplotype diplotype = geneReport.getRecommendationDiplotypes().stream().findFirst().orElse(null);
    assertNotNull(diplotype);
    assertEquals(diplotypeString, diplotype.toString());
  }

  @Test
  void testOutsideCyp2c19() {
    String gene = "CYP2C19";
    String outsideCallData = "CYP2C19\t*1/*6";
    String displayDiplotype = "*1/*6";

    OutsideCall outsideCall = new OutsideCall(s_env, outsideCallData, 0);
    GeneReport geneReport = new GeneReport(outsideCall, s_env);

    assertEquals(gene, geneReport.getGene());
    assertTrue(geneReport.isReportable());
    assertEquals(1, geneReport.getRecommendationDiplotypes().size());
    List<String> geneCalls = ReportHelpers.amdGeneCalls(geneReport);
    assertEquals(1, geneCalls.size());
    assertTrue(geneCalls.contains(displayDiplotype));
  }

  @Test
  void testNoFunctionCyp2D6() {
    OutsideCall outsideCall = new OutsideCall(s_env, "CYP2D6\t*1/*XXX", 0);

    GeneReport geneReport = new GeneReport(outsideCall, s_env);

    assertEquals("CYP2D6", geneReport.getGene());
    assertTrue(geneReport.isReportable());
    assertEquals(1, geneReport.getRecommendationDiplotypes().size());
    List<String> geneCalls = ReportHelpers.amdGeneCalls(geneReport);
    assertEquals(1, geneCalls.size());
    assertTrue(geneCalls.contains("*1/*XXX"));
    assertEquals(1, geneReport.getRecommendationDiplotypes().size());
  }

  /**
   * Regression test: {@code compareTo()} must stay consistent with {@code equals()}. Two {@link GeneReport}s
   * whose gene symbols differ only by case are not {@code equals()} (case-sensitive), so they must not
   * {@code compareTo()} as equal either.
   */
  @Test
  void testCompareToConsistentWithEqualsForDifferentCaseGene() {
    GeneReport geneReport1 = new GeneReport("CYP2D6", "test");
    GeneReport geneReport2 = new GeneReport("cyp2d6", "test");

    assertNotEquals(geneReport1, geneReport2);
    assertNotEquals(0, geneReport1.compareTo(geneReport2), "compareTo() must be consistent with equals()");

    TreeSet<GeneReport> geneReports = new TreeSet<>();
    geneReports.add(geneReport1);
    geneReports.add(geneReport2);
    assertEquals(2, geneReports.size(), "both gene reports should be retained in a TreeSet");
  }
}

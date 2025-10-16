package org.pharmgkb.pharmcat.reporter.model.result;

import java.util.List;
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
}

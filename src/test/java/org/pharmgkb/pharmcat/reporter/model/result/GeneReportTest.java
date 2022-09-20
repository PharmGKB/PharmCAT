package org.pharmgkb.pharmcat.reporter.model.result;

import java.util.SortedSet;
import java.util.TreeSet;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import org.pharmgkb.pharmcat.Env;
import org.pharmgkb.pharmcat.reporter.model.DataSource;
import org.pharmgkb.pharmcat.reporter.model.OutsideCall;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.contains;
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
    GeneReport geneReport = new GeneReport(gene, DataSource.CPIC, "test");
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

    OutsideCall outsideCall = new OutsideCall(outsideCallData);
    GeneReport geneReport = new GeneReport(outsideCall, s_env, DataSource.CPIC);

    assertEquals(gene, geneReport.getGene());
    assertTrue(geneReport.isReportable());
    assertEquals(1, geneReport.getRecommendationDiplotypes().size());
    assertEquals(1, geneReport.printDisplayCalls().size());
    assertTrue(geneReport.printDisplayCalls().contains(displayDiplotype));

    Diplotype diplotype = geneReport.getRecommendationDiplotypes().stream().findFirst().orElse(null);
    assertNotNull(diplotype);
    assertEquals(diplotypeString, diplotype.toString());
  }

  @Test
  void testOutsideCyp2c19() {
    String gene = "CYP2C19";
    String outsideCallData = "CYP2C19\t*1/*6";
    String displayDiplotype = "*1/*6";

    OutsideCall outsideCall = new OutsideCall(outsideCallData);
    GeneReport geneReport = new GeneReport(outsideCall, s_env, DataSource.CPIC);

    assertEquals(gene, geneReport.getGene());
    assertTrue(geneReport.isReportable());
    assertEquals(1, geneReport.getRecommendationDiplotypes().size());
    assertEquals(1, geneReport.printDisplayCalls().size());
    assertTrue(geneReport.printDisplayCalls().contains(displayDiplotype));
  }

  @Test
  void testNoFunctionCyp2D6() {
    OutsideCall outsideCall = new OutsideCall("CYP2D6\t*1/*XXX");

    GeneReport geneReport = new GeneReport(outsideCall, s_env, DataSource.CPIC);

    assertEquals("CYP2D6", geneReport.getGene());
    assertTrue(geneReport.isReportable());
    assertEquals(1, geneReport.getRecommendationDiplotypes().size());
    assertEquals(1, geneReport.printDisplayCalls().size());
    assertTrue(geneReport.printDisplayCalls().contains("*1/*XXX"));
    assertEquals(1, geneReport.printDisplayPhenotypes().size());
  }


  @Test
  void compare() {
    GeneReport geneReport1 = new GeneReport("TEST", DataSource.CPIC, "test");
    GeneReport geneReport2 = new GeneReport("TEST", DataSource.DPWG, "test");

    OutsideCall outsideCall = new OutsideCall("TEST\t*1/*1");
    GeneReport geneReport3 = new GeneReport(outsideCall, s_env, DataSource.CPIC);

    SortedSet<GeneReport> set = new TreeSet<>();
    set.add(geneReport1);
    set.add(geneReport2);
    set.add(geneReport3);

    assertEquals(3, set.size());
    assertThat(set, contains(geneReport3, geneReport1, geneReport2));
  }
}

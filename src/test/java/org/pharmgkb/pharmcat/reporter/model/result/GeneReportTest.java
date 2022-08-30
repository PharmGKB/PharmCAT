package org.pharmgkb.pharmcat.reporter.model.result;

import java.util.SortedSet;
import java.util.TreeSet;
import org.junit.jupiter.api.Test;
import org.pharmgkb.pharmcat.definition.PhenotypeMap;
import org.pharmgkb.pharmcat.definition.ReferenceAlleleMap;
import org.pharmgkb.pharmcat.reporter.model.DataSource;
import org.pharmgkb.pharmcat.reporter.model.OutsideCall;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.contains;
import static org.junit.jupiter.api.Assertions.*;

class GeneReportTest {
  private static final ReferenceAlleleMap sf_referenceAlleleMap = new ReferenceAlleleMap();
  private static final PhenotypeMap sf_phenotypeMap = new PhenotypeMap();

  private static final String GENE_SYMBOL1 = "GENEX";
  @Test
  void testPlaceholderGene() {
    GeneReport geneReport = new GeneReport(GENE_SYMBOL1, DataSource.CPIC, "test");
    assertEquals(GENE_SYMBOL1, geneReport.getGene());
    assertFalse(geneReport.isCalled());
    assertEquals(0, geneReport.getReporterDiplotypes().size());
  }

  private static final String GENE_SYMBOL2 = "UGT1A1";
  private static final String OUTSIDE_CALL_DATA2 = "UGT1A1\t*1/*6";
  private static final String PRINT_DIPLOTYPE2 = "*1/*6";
  private static final String DIPLOTYPE_STRING2 = "UGT1A1:*1/*6";
  @Test
  void testOutsideGene() {
    OutsideCall outsideCall = new OutsideCall(OUTSIDE_CALL_DATA2);

    GeneReport geneReport = new GeneReport(outsideCall, sf_referenceAlleleMap.get(GENE_SYMBOL2), sf_phenotypeMap,
        DataSource.CPIC);

    assertEquals(GENE_SYMBOL2, geneReport.getGene());
    assertTrue(geneReport.isReportable());
    assertEquals(1, geneReport.getReporterDiplotypes().size());
    assertEquals(1, geneReport.printDisplayCalls().size());
    assertTrue(geneReport.printDisplayCalls().contains(PRINT_DIPLOTYPE2));

    Diplotype diplotype = geneReport.getReporterDiplotypes().stream().findFirst().orElse(null);
    assertNotNull(diplotype);
    assertEquals(DIPLOTYPE_STRING2, diplotype.toString());
  }

  private static final String GENE_SYMBOL3 = "CYP2C19";
  private static final String OUTSIDE_CALL_DATA3 = "CYP2C19\t*1/*6";
  private static final String PRINT_DIPLOTYPE3 = "*1/*6";
  @Test
  void testOutsideCyp2c19() {
    OutsideCall outsideCall = new OutsideCall(OUTSIDE_CALL_DATA3);

    GeneReport geneReport = new GeneReport(outsideCall, sf_referenceAlleleMap.get(GENE_SYMBOL3), sf_phenotypeMap,
        DataSource.CPIC);

    assertEquals(GENE_SYMBOL3, geneReport.getGene());
    assertTrue(geneReport.isReportable());
    assertEquals(1, geneReport.getReporterDiplotypes().size());
    assertEquals(1, geneReport.printDisplayCalls().size());
    assertTrue(geneReport.printDisplayCalls().contains(PRINT_DIPLOTYPE3));
  }

  @Test
  void testNoFunctionCyp2D6() {
    OutsideCall outsideCall = new OutsideCall("CYP2D6\t*1/*XXX");

    GeneReport geneReport = new GeneReport(outsideCall, sf_referenceAlleleMap.get("CYP2D6"), sf_phenotypeMap,
        DataSource.CPIC);

    assertEquals("CYP2D6", geneReport.getGene());
    assertTrue(geneReport.isReportable());
    assertEquals(1, geneReport.getReporterDiplotypes().size());
    assertEquals(1, geneReport.printDisplayCalls().size());
    assertTrue(geneReport.printDisplayCalls().contains("*1/*XXX"));
    assertEquals(1, geneReport.printDisplayPhenotypes().size());
  }


  @Test
  void compare() {
    GeneReport geneReport1 = new GeneReport("TEST", DataSource.CPIC, "test");
    GeneReport geneReport2 = new GeneReport("TEST", DataSource.DPWG, "test");

    OutsideCall outsideCall = new OutsideCall("TEST\t*1/*1");
    GeneReport geneReport3 = new GeneReport(outsideCall, sf_referenceAlleleMap.get(GENE_SYMBOL3), sf_phenotypeMap,
        DataSource.CPIC);

    SortedSet<GeneReport> set = new TreeSet<>();
    set.add(geneReport1);
    set.add(geneReport2);
    set.add(geneReport3);

    assertEquals(3, set.size());
    assertThat(set, contains(geneReport3, geneReport1, geneReport2));
  }
}

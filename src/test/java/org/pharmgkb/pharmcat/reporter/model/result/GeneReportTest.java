package org.pharmgkb.pharmcat.reporter.model.result;

import org.junit.jupiter.api.Test;
import org.pharmgkb.pharmcat.definition.PhenotypeMap;
import org.pharmgkb.pharmcat.definition.ReferenceAlleleMap;
import org.pharmgkb.pharmcat.reporter.DiplotypeFactory;
import org.pharmgkb.pharmcat.reporter.model.OutsideCall;

import static org.junit.jupiter.api.Assertions.*;

class GeneReportTest {
  private static final ReferenceAlleleMap sf_referenceAlleleMap = new ReferenceAlleleMap();
  private static final PhenotypeMap sf_phenotypeMap = new PhenotypeMap();

  private static final String GENE_SYMBOL1 = "GENEX";
  @Test
  void testPlaceholderGene() {
    GeneReport geneReport = new GeneReport(GENE_SYMBOL1);
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
    DiplotypeFactory diplotypeFactory = new DiplotypeFactory(
        GENE_SYMBOL2,
        sf_phenotypeMap.lookup(GENE_SYMBOL2).orElse(null),
        sf_referenceAlleleMap.get(GENE_SYMBOL2));

    OutsideCall outsideCall = new OutsideCall(OUTSIDE_CALL_DATA2);

    GeneReport geneReport = new GeneReport(outsideCall);
    geneReport.setDiplotypes(diplotypeFactory, outsideCall);

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
    DiplotypeFactory diplotypeFactory = new DiplotypeFactory(
        GENE_SYMBOL3,
        sf_phenotypeMap.lookup(GENE_SYMBOL3).orElse(null),
        sf_referenceAlleleMap.get(GENE_SYMBOL3));

    OutsideCall outsideCall = new OutsideCall(OUTSIDE_CALL_DATA3);

    GeneReport geneReport = new GeneReport(outsideCall);
    geneReport.setDiplotypes(diplotypeFactory, outsideCall);

    assertEquals(GENE_SYMBOL3, geneReport.getGene());
    assertTrue(geneReport.isReportable());
    assertEquals(1, geneReport.getReporterDiplotypes().size());
    assertEquals(1, geneReport.printDisplayCalls().size());
    assertTrue(geneReport.printDisplayCalls().contains(PRINT_DIPLOTYPE3));
  }

  @Test
  void testNoFunctionCyp2D6() {
    DiplotypeFactory diplotypeFactory = new DiplotypeFactory(
        "CYP2D6",
        sf_phenotypeMap.lookup("CYP2D6").orElse(null),
        sf_referenceAlleleMap.get("CYP2D6"));

    OutsideCall outsideCall = new OutsideCall("CYP2D6\t*1/*XXX");

    GeneReport geneReport = new GeneReport(outsideCall);
    geneReport.setDiplotypes(diplotypeFactory, outsideCall);

    assertEquals("CYP2D6", geneReport.getGene());
    assertTrue(geneReport.isReportable());
    assertEquals(1, geneReport.getReporterDiplotypes().size());
    assertEquals(1, geneReport.printDisplayCalls().size());
    assertTrue(geneReport.printDisplayCalls().contains("*1/*XXX"));
    assertEquals(1, geneReport.printDisplayPhenotypes().size());
  }
}

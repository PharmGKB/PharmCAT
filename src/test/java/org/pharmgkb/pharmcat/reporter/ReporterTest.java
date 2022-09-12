package org.pharmgkb.pharmcat.reporter;

import org.junit.jupiter.api.Test;
import org.pharmgkb.common.util.PathUtils;
import org.pharmgkb.pharmcat.phenotype.Phenotyper;
import org.pharmgkb.pharmcat.reporter.model.DataSource;
import org.pharmgkb.pharmcat.reporter.model.result.DrugReport;
import org.pharmgkb.pharmcat.reporter.model.result.GeneReport;
import org.pharmgkb.pharmcat.reporter.model.result.GuidelineReport;

import static org.junit.jupiter.api.Assertions.*;


/**
 * Test run the Reporter and check generated data
 *
 * @author Ryan Whaley
 */
class ReporterTest {
  private static final String PHENOTYPER_FILE_PATH = "org/pharmgkb/pharmcat/phenotyper_output.json";

  @Test
  void testCypc2c9VariantPassthrough() throws Exception {

    ReportContext reportContext =
        new ReportContext(Phenotyper.read(PathUtils.getPathToResource(PHENOTYPER_FILE_PATH)).getGeneReports(), null);

    // test the CYP2C9 data
    GeneReport geneReport = reportContext.getGeneReport(DataSource.CPIC, "CYP2C9");
    assertTrue(geneReport.isReportable());
    assertTrue(geneReport.isCalled());
    assertFalse(geneReport.isOutsideCall());
    assertNotNull(geneReport.getVariantReports());
    assertTrue(
        geneReport.getVariantOfInterestReports().stream()
            .anyMatch(r -> r.getDbSnpId() != null && r.getDbSnpId().equals("rs12777823")),
        "Exemption variant not included in gene report"
    );

    // test that messages were applied for a drug
    DrugReport warfarinReport = reportContext.getDrugReports().get(DataSource.CPIC).values().stream()
        .filter(d -> d.getRelatedDrugs().contains("warfarin")).findFirst()
        .orElseThrow(() -> new RuntimeException("No warfarin drug report found"));
    assertEquals(2, warfarinReport.getMessages().size());

    // test that recommendations were matched
    DrugReport desfluraneReport = reportContext.getDrugReports().get(DataSource.CPIC).values().stream()
        .filter(d -> d.getRelatedDrugs().contains("desflurane")).findFirst()
        .orElseThrow(() -> new RuntimeException("No desflurane drug report found"));
    assertEquals(1, desfluraneReport.getGuidelines().stream().filter(GuidelineReport::isMatched).count());
  }
}

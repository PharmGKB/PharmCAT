package org.pharmgkb.pharmcat.reporter;

import org.junit.jupiter.api.Test;
import org.pharmgkb.common.util.PathUtils;
import org.pharmgkb.pharmcat.Env;
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
  private static final String PHENOTYPER_FILE_PATH = "org/pharmgkb/pharmcat/reporter/phenotyper_output.json";

  @Test
  void testCypc2c9VariantPassthrough() throws Exception {

    Phenotyper phenotyper = Phenotyper.read(PathUtils.getPathToResource(PHENOTYPER_FILE_PATH));
    ReportContext reportContext = new ReportContext(new Env(), phenotyper.getGeneReports(), null);

    // test the CYP2C9 data
    GeneReport geneReport = reportContext.getGeneReport(DataSource.CPIC, "CYP2C9");
    assertNotNull(geneReport);
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
    DrugReport warfarinReport = reportContext.getDrugReports().get(DataSource.CPIC).get("warfarin");
    assertNotNull(warfarinReport, "Missing warfarin drug report");
    assertEquals(2, warfarinReport.getMessages().size());
    assertEquals(1, warfarinReport.getGuidelines().size());
    GuidelineReport guidelineReport = warfarinReport.getGuidelines().first();
    assertEquals(1, guidelineReport.getAnnotations().size());
    assertEquals(0, guidelineReport.getAnnotations().first().getMessages().size());

    // test that recommendations were matched
    DrugReport desfluraneReport = reportContext.getDrugReports().get(DataSource.CPIC).values().stream()
        .filter(d -> d.getName().equals("desflurane")).findFirst()
        .orElseThrow(() -> new RuntimeException("No desflurane drug report found"));
    assertEquals(1, desfluraneReport.getGuidelines().stream().filter(GuidelineReport::isMatched).count());
  }
}

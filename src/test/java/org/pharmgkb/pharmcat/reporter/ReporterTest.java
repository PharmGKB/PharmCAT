package org.pharmgkb.pharmcat.reporter;

import java.nio.file.Files;
import java.nio.file.Path;
import org.jsoup.Jsoup;
import org.jsoup.nodes.Document;
import org.jsoup.select.Elements;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.TestInfo;
import org.pharmgkb.common.util.PathUtils;
import org.pharmgkb.pharmcat.BaseConfig;
import org.pharmgkb.pharmcat.Env;
import org.pharmgkb.pharmcat.TestUtils;
import org.pharmgkb.pharmcat.phenotype.Phenotyper;
import org.pharmgkb.pharmcat.reporter.format.HtmlFormat;
import org.pharmgkb.pharmcat.reporter.handlebars.ReportHelpers;
import org.pharmgkb.pharmcat.reporter.model.DataSource;
import org.pharmgkb.pharmcat.reporter.model.PrescribingGuidanceSource;
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
  private Env m_env;

  @BeforeEach
  void before() throws Exception {
    ReportHelpers.setDebugMode(true);
    m_env = new Env();
    //TestUtils.setSaveTestOutput(true);
  }

  @AfterEach
  void deleteDirectory(TestInfo testInfo) {
    TestUtils.deleteTestOutputDirectory(testInfo);
  }


  @Test
  void cypc2c9VariantPassthrough() throws Exception {

    Phenotyper phenotyper = Phenotyper.read(PathUtils.getPathToResource("org/pharmgkb/pharmcat/reporter/ReporterTest-cypc2c9VariantPassthrough.json"));
    ReportContext reportContext = new ReportContext(m_env, phenotyper.getGeneReports(), null);

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
    DrugReport warfarinReport = reportContext.getDrugReports().get(PrescribingGuidanceSource.CPIC_GUIDELINE).get("warfarin");
    assertNotNull(warfarinReport, "Missing warfarin drug report");
    assertEquals(2, warfarinReport.getMessages().size());
    assertEquals(1, warfarinReport.getGuidelines().size());
    GuidelineReport guidelineReport = warfarinReport.getGuidelines().first();
    assertEquals(1, guidelineReport.getAnnotations().size());
    assertEquals(0, guidelineReport.getAnnotations().first().getMessages().size());

    // test that recommendations were matched
    DrugReport desfluraneReport = reportContext.getDrugReport(PrescribingGuidanceSource.CPIC_GUIDELINE, "desflurane");
    assertNotNull(desfluraneReport);
    assertEquals(1, desfluraneReport.getGuidelines().stream().filter(GuidelineReport::isMatched).count());
  }

  @Test
  void multipleActivityScores(TestInfo testInfo) throws Exception {

    Phenotyper phenotyper = Phenotyper.read(
        PathUtils.getPathToResource("org/pharmgkb/pharmcat/reporter/ReporterTest-multipleActivityScores.json"));
    ReportContext reportContext = new ReportContext(new Env(), phenotyper.getGeneReports(), null);

    // test the CYP2C9 data
    GeneReport geneReport = reportContext.getGeneReport(DataSource.CPIC, "CYP2C9");
    assertNotNull(geneReport);
    assertTrue(geneReport.isReportable());
    assertTrue(geneReport.isCalled());
    assertFalse(geneReport.isOutsideCall());
    assertNotNull(geneReport.getVariantReports());

    // CYP2C9 has 2 activity scores - should have 2 rows (AnnotationReports)
    DrugReport tenoxicam = reportContext.getDrugReport(PrescribingGuidanceSource.CPIC_GUIDELINE, "tenoxicam");
    assertNotNull(tenoxicam);
    assertEquals(1, tenoxicam.getGuidelines().size());
    assertEquals(2, tenoxicam.getGuidelines().first().getAnnotations().size());
    assertEquals(1, tenoxicam.getGuidelines().first().getAnnotations().first().getActivityScores().size());
    assertEquals(1, tenoxicam.getGuidelines().first().getAnnotations().first().getPhenotypes().size());

    // warfarin is a special case - even though CYP2C9 has 2 activity scores, it gets merged into 1 row/AnnotationReport
    DrugReport warfarin = reportContext.getDrugReport(PrescribingGuidanceSource.CPIC_GUIDELINE, "warfarin");
    assertNotNull(warfarin);
    assertEquals(1, warfarin.getGuidelines().size());
    assertEquals(1, warfarin.getGuidelines().first().getAnnotations().size());
    assertEquals(0, warfarin.getGuidelines().first().getAnnotations().first().getActivityScores().size());
    assertEquals(0, warfarin.getGuidelines().first().getAnnotations().first().getPhenotypes().size());

    Path reporterOutput = printReport(testInfo, reportContext);
    Document document = Jsoup.parse(reporterOutput.toFile());
    // should have no tags for CPIC
    Elements tags = document.select(".cpic-guideline-warfarin .tag");
    assertEquals(0, tags.size());
    // should have tags for DPWG
    tags = document.select(".dpwg-guideline-warfarin .tag");
    assertEquals(2, tags.size());

    Elements phenotype = document.select(".cpic-guideline-warfarin .rx-phenotype");
    assertEquals(0, phenotype.size());
    phenotype = document.select(".dpwg-guideline-warfarin .rx-phenotype");
    assertEquals(2, phenotype.size());

    Elements activityScore = document.select(".cpic-guideline-warfarin .rx-activity");
    assertEquals(0, activityScore.size());
    activityScore = document.select(".dpwg-guideline-warfarin .rx-activity");
    assertEquals(1, activityScore.size());
  }

  @Test
  void multiplePhenotypes(TestInfo testInfo) throws Exception {

    Phenotyper phenotyper = Phenotyper.read(
        PathUtils.getPathToResource("org/pharmgkb/pharmcat/reporter/ReporterTest-multiplePhenotypes.json"));
    ReportContext reportContext = new ReportContext(new Env(), phenotyper.getGeneReports(), null);

    // test the CYP2C9 data
    GeneReport geneReport = reportContext.getGeneReport(DataSource.CPIC, "CYP2C9");
    assertNotNull(geneReport);
    assertTrue(geneReport.isReportable());
    assertTrue(geneReport.isCalled());
    assertFalse(geneReport.isOutsideCall());
    assertNotNull(geneReport.getVariantReports());

    // CYP2C9 has 2 activity scores - should have 2 rows (AnnotationReports)
    DrugReport tenoxicam = reportContext.getDrugReport(PrescribingGuidanceSource.CPIC_GUIDELINE, "tenoxicam");
    assertNotNull(tenoxicam);
    assertEquals(1, tenoxicam.getGuidelines().size());
    assertEquals(2, tenoxicam.getGuidelines().first().getAnnotations().size());
    assertEquals(1, tenoxicam.getGuidelines().first().getAnnotations().first().getActivityScores().size());
    assertEquals(1, tenoxicam.getGuidelines().first().getAnnotations().first().getPhenotypes().size());

    // warfarin is a special case - even though CYP2C9 has 2 phenotypes, it gets merged into 1 row/AnnotationReport
    DrugReport warfarin = reportContext.getDrugReport(PrescribingGuidanceSource.CPIC_GUIDELINE, "warfarin");
    assertNotNull(warfarin);
    assertEquals(1, warfarin.getGuidelines().size());
    assertEquals(1, warfarin.getGuidelines().first().getAnnotations().size());
    assertEquals(0, warfarin.getGuidelines().first().getAnnotations().first().getActivityScores().size());
    assertEquals(0, warfarin.getGuidelines().first().getAnnotations().first().getPhenotypes().size());

    Path reporterOutput = printReport(testInfo, reportContext);
    Document document = Jsoup.parse(reporterOutput.toFile());
    // should have no tags for CPIC
    Elements tags = document.select(".cpic-guideline-warfarin .tag");
    assertEquals(0, tags.size());
    // should have tags for DPWG
    tags = document.select(".dpwg-guideline-warfarin .tag");
    assertEquals(1, tags.size());

    Elements phenotype = document.select(".cpic-guideline-warfarin .rx-phenotype");
    assertEquals(0, phenotype.size());
    phenotype = document.select(".dpwg-guideline-warfarin .rx-phenotype");
    assertEquals(1, phenotype.size());

    Elements activityScore = document.select(".cpic-guideline-warfarin .rx-activity");
    assertEquals(0, activityScore.size());
    activityScore = document.select(".dpwg-guideline-warfarin .rx-activity");
    assertEquals(1, activityScore.size());
  }

  private Path printReport(TestInfo testInfo, ReportContext reportContext) throws Exception {

    Path outputDir = TestUtils.getTestOutputDir(testInfo, false);
    if (!Files.isDirectory(outputDir)) {
      Files.createDirectories(outputDir);
    }
    Path file = outputDir.resolve("test" + BaseConfig.REPORTER_SUFFIX + ".html");
    if (TestUtils.isSaveTestOutput()) {
      System.out.println("Saving report to " + file);
    }
    new HtmlFormat(file, m_env, true)
        .compact(true)
        .write(reportContext);
    return file;
  }
}

package org.pharmgkb.pharmcat;

import java.util.List;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.TestInfo;
import org.pharmgkb.pharmcat.reporter.handlebars.ReportHelpers;
import org.pharmgkb.pharmcat.reporter.model.PrescribingGuidanceSource;
import org.pharmgkb.pharmcat.reporter.model.result.DrugReport;
import org.pharmgkb.pharmcat.reporter.model.result.Genotype;
import org.pharmgkb.pharmcat.reporter.model.result.GuidelineReport;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;


/**
 * These are JUnit tests for toxic genes.
 *
 * @author Mark Woon
 */
class ToxicGenesTest {


  @BeforeAll
  static void prepare() {
    ReportHelpers.setDebugMode(true);
    TestUtils.setSaveTestOutput(true);
  }

  @AfterEach
  void deleteDirectory(TestInfo testInfo) {
    TestUtils.deleteTestOutputDirectory(testInfo);
  }


  @Test
  void testNudt15Ref(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("NUDT15");
    testWrapper.execute();

    testWrapper.testPrintCpicCalls("NUDT15", "*1/*1");
    testWrapper.testRecommendedDiplotypes("NUDT15", "*1", "*1");

    testWrapper.testMatchedAnnotations("azathioprine", 2);
    testWrapper.testMatchedAnnotations("mercaptopurine", 2);
    testWrapper.testAnyMatchFromSource("mercaptopurine", PrescribingGuidanceSource.CPIC_GUIDELINE);
    testWrapper.testAnyMatchFromSource("mercaptopurine", PrescribingGuidanceSource.DPWG_GUIDELINE);
    testWrapper.testMatchedAnnotations("thioguanine", 2);
    testWrapper.testAnyMatchFromSource("thioguanine", PrescribingGuidanceSource.CPIC_GUIDELINE);
    testWrapper.testAnyMatchFromSource("thioguanine", PrescribingGuidanceSource.DPWG_GUIDELINE);
  }

  @Test
  void testNudt15S2(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("NUDT15", "rs746071566", "GAGTCG(3)", "GAGTCG(4)")
        .variation("NUDT15", "rs116855232", "C", "T")
    ;
    testWrapper.execute();

    testWrapper.testCalledByMatcher("NUDT15");
    testWrapper.testPrintCpicCalls("NUDT15", "*1/*2");
    testWrapper.testRecommendedDiplotypes("NUDT15", "*1", "*2");

    DrugReport azaReport = testWrapper.getContext().getDrugReport(PrescribingGuidanceSource.CPIC_GUIDELINE, "azathioprine");
    assertNotNull(azaReport);
    GuidelineReport azaCpicGuideline = azaReport.getGuidelines().iterator().next();
    List<Genotype> genotypes = Genotype.makeGenotypes(azaCpicGuideline.getGeneReports());
    assertEquals(1, genotypes.size());

    testWrapper.testMatchedAnnotations("azathioprine", 4);
  }

  @Test
  void testNudt15S3(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("NUDT15", "rs116855232", "C", "T")
    ;
    testWrapper.execute();

    testWrapper.testCalledByMatcher("NUDT15");
    testWrapper.testPrintCpicCalls("NUDT15", "*1/*3");
    testWrapper.testRecommendedDiplotypes("NUDT15", "*1", "*3");

    testWrapper.testMatchedAnnotations("azathioprine", 4);
    testWrapper.testMatchedAnnotations("mercaptopurine", 4);
    testWrapper.testMatchedAnnotations("thioguanine", 4);
  }
}

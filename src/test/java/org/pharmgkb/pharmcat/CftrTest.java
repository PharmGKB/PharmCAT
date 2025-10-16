package org.pharmgkb.pharmcat;

import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import org.jsoup.nodes.Document;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.TestInfo;
import org.pharmgkb.pharmcat.reporter.TextConstants;
import org.pharmgkb.pharmcat.reporter.format.html.ReportHelpers;

import static org.pharmgkb.pharmcat.PipelineTest.*;


/**
 * This is a JUnit test for CFTR.
 *
 * @author Mark Woon
 */
class CftrTest {

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
  void ref_ref(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("CFTR");
    Path vcfFile = testWrapper.execute();

    List<String> expectedCalls = List.of(TextConstants.HOMOZYGOUS_REFERENCE);

    testWrapper.testCalledByMatcher("CFTR");
    testWrapper.testSourceDiplotypes("CFTR", expectedCalls);
    testWrapper.testRecommendedDiplotypes("CFTR",
        List.of("ivacaftor non-responsive CFTR sequence", "ivacaftor non-responsive CFTR sequence"));
    testWrapper.testPrintCalls("CFTR", expectedCalls);

    testWrapper.testMatchedAnnotations("ivacaftor", 1);

    Document document = readHtmlReport(vcfFile);
    htmlChecks(document, "CFTR", expectedCalls, "ivacaftor", RecPresence.YES, RecPresence.NO);
  }

  @Test
  void d1270NHet(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("CFTR", "rs11971167", "G", "A");
    Path vcfFile = testWrapper.execute();

    List<String> expectedCalls = List.of("D1270N (heterozygous)");

    testWrapper.testCalledByMatcher("CFTR");
    testWrapper.testSourceDiplotypes("CFTR", expectedCalls);
    testWrapper.testRecommendedDiplotypes("CFTR", List.of("ivacaftor non-responsive CFTR sequence", "D1270N"));
    testWrapper.testPrintCalls("CFTR", expectedCalls);

    testWrapper.testMatchedAnnotations("ivacaftor", 2);

    Document document = readHtmlReport(vcfFile);
    htmlChecks(document, "CFTR", expectedCalls, "ivacaftor", RecPresence.YES, RecPresence.NO);
  }

  @Test
  void d1270NG551D(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("CFTR", "rs11971167", "G", "A")
        .variation("CFTR", "rs75527207", "G", "A");
    Path vcfFile = testWrapper.execute();

    List<String> expectedCalls = List.of("D1270N/G551D");

    testWrapper.testCalledByMatcher("CFTR");
    testWrapper.testSourceDiplotypes("CFTR", expectedCalls);
    testWrapper.testRecommendedDiplotypes("CFTR", expectedCallsToRecommendedDiplotypes(expectedCalls));
    testWrapper.testPrintCalls("CFTR", expectedCalls);

    testWrapper.testMatchedAnnotations("ivacaftor", 2);

    Document document = readHtmlReport(vcfFile);
    htmlChecks(document, "CFTR", expectedCalls, "ivacaftor", RecPresence.YES, RecPresence.NO);
  }


  /**
   *  Check to see if the CPIC-specified allele name of "ivacaftor non-responsive CFTR sequence" gets translated to the
   *  more generic value of "Reference" for display in the report
   */
  @Test
  void outsideReferenceCall(TestInfo testInfo) throws Exception {
    Path outsideCallPath = TestUtils.createTestFile(testInfo, ".tsv");
    try (PrintWriter writer = new PrintWriter(Files.newBufferedWriter(outsideCallPath))) {
      writer.println("CFTR\tivacaftor non-responsive CFTR sequence/ivacaftor non-responsive CFTR sequence");
    }
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("CYP2C19")
        .reference("CYP2C9")
    ;
    testWrapper.executeWithOutsideCalls(outsideCallPath);

    testWrapper.testCalledByMatcher("CYP2C19");
    testWrapper.testCalledByMatcher("CYP2C9");

    testWrapper.testReportable("CFTR");
    testWrapper.testPrintCalls("CFTR", "Reference/Reference");
  }
}

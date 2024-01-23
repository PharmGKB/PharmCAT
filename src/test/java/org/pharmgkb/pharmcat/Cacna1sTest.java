package org.pharmgkb.pharmcat;

import java.nio.file.Path;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import org.jsoup.nodes.Document;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.TestInfo;
import org.pharmgkb.pharmcat.reporter.TextConstants;
import org.pharmgkb.pharmcat.reporter.handlebars.ReportHelpers;
import org.pharmgkb.pharmcat.reporter.model.DataSource;

import static org.pharmgkb.pharmcat.PipelineTest.*;


/**
 * This JUnit test validates CACNA1S through the pipeline.
 *
 * @author Mark Woon
 */
class Cacna1sTest {

  @BeforeAll
  static void prepare() {
    ReportHelpers.setDebugMode(true);
    //TestUtils.setSaveTestOutput(true);
  }

  @AfterEach
  void deleteDirectory(TestInfo testInfo) {
    TestUtils.deleteTestOutputDirectory(testInfo);
  }


  /**
   * Test homozygous CACNA1S with 2 het RYR1 calls.
   */
  @Test
  void cacna1sHom_ryr1Het2(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false, false, false);
    testWrapper.getVcfBuilder()
        .reference("CACNA1S")
        .variation("RYR1", "rs193922746", "A", "G") // c.97A>G
        .variation("RYR1", "rs193922749", "C", "A") // c.152C>A
    ;
    Path vcfFile = testWrapper.execute(null);

    testWrapper.testCalledByMatcher("CACNA1S");
    testWrapper.testCalledByMatcher("RYR1");

    List<String> cacna1sExpectedCalls = List.of(TextConstants.REFERENCE + TextConstants.GENOTYPE_DELIMITER + TextConstants.REFERENCE);
    testWrapper.testSourceDiplotypes(DataSource.CPIC, "CACNA1S", cacna1sExpectedCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "CACNA1S", expectedCallsToRecommendedDiplotypes(cacna1sExpectedCalls));

    List<String> ryr1ExpectedCalls = List.of("c.97A>G", "c.152C>A");
    testWrapper.testSourceDiplotypes(DataSource.CPIC, "RYR1", ryr1ExpectedCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "RYR1", ryr1ExpectedCalls);

    Document document = readHtmlReport(vcfFile);
    Map<String, List<String>> expectedCallsMap = new LinkedHashMap<>();
    expectedCallsMap.put("CACNA1S", cacna1sExpectedCalls);
    expectedCallsMap.put("RYR1", ryr1ExpectedCalls);
    Map<String, List<String>> cpicStyleCallsMap = new LinkedHashMap<>();
    cpicStyleCallsMap.put("RYR1", ryr1ExpectedCalls);
    htmlChecks(document, expectedCallsMap, cpicStyleCallsMap, "enflurane", RecPresence.YES, RecPresence.NO);
  }


  /**
   * Test a het CACNA1S and missing RYR1 call
   */
  @Test
  void cacna1sHet_ryr1Missing(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false, false, false);
    testWrapper.getVcfBuilder()
        .reference("CACNA1S")
        .variation("CACNA1S", "rs772226819", "G", "A");
    Path vcfFile = testWrapper.execute(null);

    testWrapper.testNotCalledByMatcher("RYR1");
    testWrapper.testCalledByMatcher("CACNA1S");

    List<String> cacna1sExpectedCalls = List.of(TextConstants.REFERENCE + TextConstants.GENOTYPE_DELIMITER + "c.520C>T");
    List<String> cacna1sCpicStyleCalls = List.of("c.520C>T (heterozygous)");

    testWrapper.testSourceDiplotypes(DataSource.CPIC, "CACNA1S", cacna1sCpicStyleCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "CACNA1S", expectedCallsToRecommendedDiplotypes(cacna1sExpectedCalls));
    testWrapper.testMatchedAnnotations("enflurane", DataSource.CPIC, 1);

    Document document = readHtmlReport(vcfFile);
    Map<String, List<String>> expectedCallsMap = new LinkedHashMap<>();
    expectedCallsMap.put("CACNA1S", cacna1sExpectedCalls);
    expectedCallsMap.put("RYR1", NO_DATA);
    Map<String, List<String>> cpicStyleCallsMap = new LinkedHashMap<>();
    cpicStyleCallsMap.put("CACNA1S", cacna1sCpicStyleCalls);
    htmlChecks(document, expectedCallsMap, cpicStyleCallsMap, "enflurane", RecPresence.YES, RecPresence.NO);
  }

  /**
   * Test a het CACNA1S and het RYR1 call.
   */
  @Test
  void cacna1sHet_ryr1Het(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false, false, false);
    testWrapper.getVcfBuilder()
        .reference("CACNA1S")
        .reference("RYR1")
        .variation("CACNA1S", "rs772226819", "G", "A")
        .variation("RYR1", "rs118192178", "G", "C");
    Path vcfFile = testWrapper.execute(null);

    testWrapper.testCalledByMatcher("CACNA1S");
    testWrapper.testCalledByMatcher("RYR1");

    List<String> ryr1ExpectedCalls = List.of("Reference/c.7522C>G");
    List<String> ryr1CpicStyleCalls = List.of("c.7522C>G (heterozygous)");

    testWrapper.testSourceDiplotypes(DataSource.CPIC, "CACNA1S", List.of("c.520C>T (heterozygous)"));
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "CACNA1S", List.of(TextConstants.REFERENCE, "c.520C>T"));

    testWrapper.testSourceDiplotypes(DataSource.CPIC, "RYR1", ryr1CpicStyleCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "RYR1", List.of(TextConstants.REFERENCE, "c.7522C>G"));

    Document document = readHtmlReport(vcfFile);
    htmlChecks(document, "CACNA1S", List.of("c.520C>T (heterozygous)"), null, RecPresence.YES, RecPresence.NO);
    htmlChecks(document, "RYR1", ryr1ExpectedCalls, ryr1CpicStyleCalls, null, RecPresence.YES, RecPresence.NO);

    Map<String, List<String>> expectedCallsMap = new LinkedHashMap<>();
    expectedCallsMap.put("CACNA1S", List.of("c.520C>T (heterozygous)"));
    expectedCallsMap.put("RYR1", ryr1ExpectedCalls);
    Map<String, List<String>> cpicStyleCallsMap = new LinkedHashMap<>();
    cpicStyleCallsMap.put("RYR1", ryr1CpicStyleCalls);
    htmlChecks(document, expectedCallsMap, cpicStyleCallsMap, "enflurane", RecPresence.YES, RecPresence.NO);
  }

  /**
   * Test a het CACNA1S with ref RYR1.
   */
  @Test
  void cacna1sHet_ryr1Ref(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false, false, false);
    testWrapper.getVcfBuilder()
        .variation("CACNA1S", "rs772226819", "G", "A")
        .reference("RYR1")
    ;
    Path vcfFile = testWrapper.execute(null);

    testWrapper.testCalledByMatcher("CACNA1S");
    testWrapper.testCalledByMatcher("RYR1");

    List<String> cacna1sExpectedCalls = List.of(TextConstants.REFERENCE, "c.520C>T");
    List<String> cacna1sCpicStyleCalls = List.of("c.520C>T (heterozygous)");
    testWrapper.testSourceDiplotypes(DataSource.CPIC, "CACNA1S", cacna1sCpicStyleCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "CACNA1S", cacna1sExpectedCalls);

    List<String> ryr1ExpectedCalls = List.of(TextConstants.REFERENCE + TextConstants.GENOTYPE_DELIMITER + TextConstants.REFERENCE);
    testWrapper.testSourceDiplotypes(DataSource.CPIC, "RYR1", ryr1ExpectedCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "RYR1", expectedCallsToRecommendedDiplotypes(ryr1ExpectedCalls));

    Document document = readHtmlReport(vcfFile);
    Map<String, List<String>> expectedCallsMap = new LinkedHashMap<>();
    expectedCallsMap.put("CACNA1S", cacna1sExpectedCalls);
    expectedCallsMap.put("RYR1", ryr1ExpectedCalls);
    Map<String, List<String>> cpicStyleCallsMap = new LinkedHashMap<>();
    cpicStyleCallsMap.put("CACNA1S", cacna1sCpicStyleCalls);
    htmlChecks(document, expectedCallsMap, cpicStyleCallsMap, "enflurane", RecPresence.YES, RecPresence.NO);
  }
}

package org.pharmgkb.pharmcat;

import java.nio.file.Path;
import java.util.List;
import java.util.SortedMap;
import java.util.TreeMap;
import com.google.common.collect.ImmutableSortedMap;
import com.google.common.collect.Ordering;
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
 * This JUnit test validates RYR1 through the pipeline.
 *
 * @author Mark Woon
 */
class Ryr1Test {

  @BeforeAll
  static void prepare() {
    ReportHelpers.setDebugMode(true);
    TestUtils.setSaveTestOutput(true);
  }

  @AfterEach
  void deleteDirectory(TestInfo testInfo) {
    TestUtils.deleteTestOutputDirectory(testInfo);
  }


  /**
   * Test a ref CACNA1S with 1 undocumented het RYR1 call.
   */
  @Test
  void cacna1sRef_ryr1Undocumented(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false, false, false);
    testWrapper.getVcfBuilder()
        .allowUnknownAllele()
        .reference("CACNA1S")
        .reference("RYR1")
        .variation("RYR1", "rs193922749", "A", "T"); // c.152C>A
    Path vcfFile = testWrapper.execute(null);

    testWrapper.testCalledByMatcher("CACNA1S");
    testWrapper.testCalledByMatcher("RYR1");

    testWrapper.testSourceDiplotypes(DataSource.CPIC, "CACNA1S", List.of(TextConstants.HOMOZYGOUS_REFERENCE));
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "CACNA1S", List.of(TextConstants.REFERENCE, TextConstants.REFERENCE));

    List<String> ryr1ExpectedCalls = List.of("Reference/c.152C>A");
    List<String> ryr1CpicStyleCalls = List.of("c.152C>A (heterozygous)");
    testWrapper.testSourceDiplotypes(DataSource.CPIC, "RYR1", ryr1ExpectedCalls, ryr1CpicStyleCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "RYR1", expectedCallsToRecommendedDiplotypes(ryr1ExpectedCalls));

    Document document = readHtmlReport(vcfFile);
    SortedMap<String, List<String>> expectedCallsMap = new TreeMap<>();
    expectedCallsMap.put("CACNA1S", List.of(TextConstants.HOMOZYGOUS_REFERENCE));
    expectedCallsMap.put("RYR1", ryr1ExpectedCalls);
    SortedMap<String, List<String>> cpicStyleCallsMap = new TreeMap<>();
    cpicStyleCallsMap.put("RYR1", ryr1CpicStyleCalls);
    htmlChecks(document, expectedCallsMap, cpicStyleCallsMap, "enflurane", RecPresence.YES, RecPresence.NO);
  }

  /**
   * Test a ref CACNA1S with 1 het RYR1 call.
   */
  @Test
  void cacna1sRef_ryr1Het(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false, false, false);
    testWrapper.getVcfBuilder()
        .reference("CACNA1S")
        .reference("RYR1")
        .variation("RYR1", "rs193922749", "A", "C"); // c.152C>A
    Path vcfFile = testWrapper.execute(null);

    testWrapper.testCalledByMatcher("CACNA1S");
    testWrapper.testCalledByMatcher("RYR1");

    testWrapper.testSourceDiplotypes(DataSource.CPIC, "CACNA1S", List.of(TextConstants.HOMOZYGOUS_REFERENCE));
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "CACNA1S", List.of(TextConstants.REFERENCE, TextConstants.REFERENCE));

    List<String> ryr1ExpectedCalls = List.of("Reference/c.152C>A");
    List<String> ryr1CpicStyleCalls = List.of("c.152C>A (heterozygous)");
    testWrapper.testSourceDiplotypes(DataSource.CPIC, "RYR1", ryr1ExpectedCalls, ryr1CpicStyleCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "RYR1", expectedCallsToRecommendedDiplotypes(ryr1ExpectedCalls));

    Document document = readHtmlReport(vcfFile);
    SortedMap<String, List<String>> expectedCallsMap = new TreeMap<>();
    expectedCallsMap.put("CACNA1S", List.of(TextConstants.HOMOZYGOUS_REFERENCE));
    expectedCallsMap.put("RYR1", ryr1ExpectedCalls);
    SortedMap<String, List<String>> cpicStyleCallsMap = new TreeMap<>();
    cpicStyleCallsMap.put("RYR1", ryr1CpicStyleCalls);
    htmlChecks(document, expectedCallsMap, cpicStyleCallsMap, "enflurane", RecPresence.YES, RecPresence.NO);
  }

  /**
   * Test a ref CACNA1S with 2 het RYR1 calls.
   */
  @Test
  void cacna1sRef_ryr1Het2(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false, false, false);
    testWrapper.getVcfBuilder()
        .reference("CACNA1S")
        .reference("RYR1")
        .variation("RYR1", "rs34694816", "A", "G")
        .variation("RYR1", "rs137933390", "A", "G")
    ;
    Path vcfFile = testWrapper.execute(null);

    testWrapper.testCalledByMatcher("CACNA1S");
    testWrapper.testCalledByMatcher("RYR1");

    List<String> ryr1ExpectedCalls = List.of("c.4024A>G", "c.4178A>G");
    testWrapper.testSourceDiplotypes(DataSource.CPIC, "RYR1", ryr1ExpectedCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "RYR1", ryr1ExpectedCalls);

    Document document = readHtmlReport(vcfFile);
    SortedMap<String, List<String>> expectedCallsMap = new TreeMap<>();
    expectedCallsMap.put("CACNA1S", List.of(TextConstants.HOMOZYGOUS_REFERENCE));
    expectedCallsMap.put("RYR1", ryr1ExpectedCalls);
    htmlChecks(document, expectedCallsMap, "enflurane", RecPresence.YES, RecPresence.NO);
  }

  /**
   * Test a ref CACNA1S with 3 het RYR1 calls.
   */
  @Test
  void cacna1sRef_ryr1Het3(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false, false, false);
    testWrapper.getVcfBuilder()
        .reference("CACNA1S")
        .reference("RYR1")
        .variation("RYR1", "rs34694816", "A", "G")   // c.4024A>G, normal
        .variation("RYR1", "rs137933390", "A", "G")  // c.4178A>G, normal
        .variation("RYR1", "rs145573319", "A", "G")  // c.4400A>G, uncertain
    ;
    Path vcfFile = testWrapper.execute(null);

    testWrapper.testCalledByMatcher("CACNA1S");
    testWrapper.testCalledByMatcher("RYR1");

    List<String> ryr1ExpectedCalls = List.of("c.4024A>G", "c.4178A>G", "c.4400A>G");
    testWrapper.testSourceDiplotypes(DataSource.CPIC, "RYR1", ryr1ExpectedCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "RYR1", List.of("c.4024A>G", "c.4178A>G"));

    Document document = readHtmlReport(vcfFile);
    SortedMap<String, List<String>> expectedCallsMap = new TreeMap<>();
    expectedCallsMap.put("CACNA1S", List.of(TextConstants.HOMOZYGOUS_REFERENCE));
    expectedCallsMap.put("RYR1", ryr1ExpectedCalls);
    htmlChecks(document, expectedCallsMap, "enflurane", RecPresence.YES, RecPresence.NO);
  }


  @Test
  void cacna1sMissing_ryr1Ref(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("RYR1");
    Path vcfFile = testWrapper.execute(null);

    testWrapper.testCalledByMatcher("RYR1");

    List<String> expectedCalls = List.of(TextConstants.HOMOZYGOUS_REFERENCE);
    testWrapper.testSourceDiplotypes(DataSource.CPIC, "RYR1", expectedCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "RYR1", List.of(TextConstants.REFERENCE, TextConstants.REFERENCE));
    testWrapper.testPrintCalls(DataSource.CPIC, "RYR1", expectedCalls);

    testWrapper.testMatchedAnnotations("desflurane", DataSource.CPIC, 1);
    testWrapper.testNoMatchFromSource("desflurane", DataSource.DPWG);

    Document document = readHtmlReport(vcfFile);
    htmlChecks(document,
        new ImmutableSortedMap.Builder<String, List<String>>(Ordering.natural())
            .put("CACNA1S", NO_DATA)
            .put("RYR1", expectedCalls)
            .build(),
        null, "desflurane", RecPresence.YES,
        new ImmutableSortedMap.Builder<String, String>(Ordering.natural())
            .put("CACNA1S", "No Result")
            .put("RYR1", "Uncertain Susceptibility")
            .build(),
        RecPresence.NO, null);
  }

  /**
   * Test a missing CACNA1S gene with 1 het RYR1 call.
   */
  @Test
  void cacna1sMissing_ryr1Het(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false, false, false);
    testWrapper.getVcfBuilder()
        .reference("RYR1")
        .variation("RYR1", "rs118192178", "C", "T");
    Path vcfFile = testWrapper.execute(null);

    testWrapper.testNotCalledByMatcher("CACNA1S");
    testWrapper.testCalledByMatcher("RYR1");

    List<String> ryr1ExpectedCalls = List.of(TextConstants.REFERENCE + TextConstants.GENOTYPE_DELIMITER + "c.7522C>T");
    List<String> ryr1CpicStyleCalls = List.of("c.7522C>T (heterozygous)");

    testWrapper.testSourceDiplotypes(DataSource.CPIC, "RYR1", ryr1ExpectedCalls, ryr1CpicStyleCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "RYR1", expectedCallsToRecommendedDiplotypes(ryr1ExpectedCalls));

    Document document = readHtmlReport(vcfFile);

    SortedMap<String, List<String>> expectedCallsMap = new TreeMap<>();
    expectedCallsMap.put("CACNA1S", NO_DATA);
    expectedCallsMap.put("RYR1", ryr1ExpectedCalls);
    SortedMap<String, List<String>> cpicStyleCallsMap = new TreeMap<>();
    cpicStyleCallsMap.put("RYR1", ryr1CpicStyleCalls);
    htmlChecks(document, expectedCallsMap, cpicStyleCallsMap, "enflurane", RecPresence.YES, RecPresence.NO);
  }

  /**
   * Test a missing CACNA1S gene with 2 RYR1 het calls.
   */
  @Test
  void cacna1sMissing_ryr1Het2(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false, false, false);
    testWrapper.getVcfBuilder()
        .reference("RYR1")
        .variation("RYR1", "rs34694816", "A", "G")
        .variation("RYR1", "rs137933390", "A", "G")  // c.4178A>G
    ;
    Path vcfFile = testWrapper.execute(null);

    testWrapper.testNotCalledByMatcher("CACNA1S");
    testWrapper.testCalledByMatcher("RYR1");

    List<String> ryr1ExpectedCalls = List.of("c.4024A>G", "c.4178A>G");
    testWrapper.testSourceDiplotypes(DataSource.CPIC, "RYR1", ryr1ExpectedCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "RYR1", ryr1ExpectedCalls);

    Document document = readHtmlReport(vcfFile);
    htmlChecks(document,
        new ImmutableSortedMap.Builder<String, List<String>>(Ordering.natural())
            .put("CACNA1S", NO_DATA)
            .put("RYR1", ryr1ExpectedCalls)
            .build(),
        null, "enflurane", RecPresence.YES,
        new ImmutableSortedMap.Builder<String, String>(Ordering.natural())
            .put("CACNA1S", "No Result")
            .put("RYR1", "Uncertain Susceptibility")
            .build(),
            RecPresence.NO, null);
  }

  /**
   * Test a missing CACNA1S gene with 2 RYR1 het calls and 1 hom call.
   */
  @Test
  void cacna1sMissing_ryr1Het2Hom1(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false, false, false);
    testWrapper.getVcfBuilder()
        .reference("RYR1")
        .variation("RYR1", "rs193922746", "G", "G") // c.97A>G
        .variation("RYR1", "rs193922749", "C", "A") // c.152C>A
        .variation("RYR1", "rs142474192", "G", "A")  // c.418G>A
    ;
    Path vcfFile = testWrapper.execute(null);

    testWrapper.testNotCalledByMatcher("CACNA1S");
    testWrapper.testCalledByMatcher("RYR1");

    List<String> ryr1ExpectedCalls = List.of("c.97A>G (homozygous)", "c.152C>A", "c.418G>A");
    testWrapper.testSourceDiplotypes(DataSource.CPIC, "RYR1", ryr1ExpectedCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "RYR1", List.of("c.97A>G", "c.97A>G"));

    Document document = readHtmlReport(vcfFile);
    SortedMap<String, List<String>> expectedCallsMap = new TreeMap<>();
    expectedCallsMap.put("CACNA1S", NO_DATA);
    expectedCallsMap.put("RYR1", ryr1ExpectedCalls);
    htmlChecks(document, expectedCallsMap, "enflurane", RecPresence.YES, RecPresence.NO);
  }

  /**
   * Test a missing CACNA1S gene with 2 RYR1 het calls and 1 hom call, phased.
   */
  @Test
  void cacna1sMissing_ryr1Het2Hom1Phased(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false, false, false);
    testWrapper.getVcfBuilder()
        .phased()
        .reference("RYR1")
        .variation("RYR1", "rs193922746", "G", "G") // c.97A>G
        .variation("RYR1", "rs193922749", "C", "A") // c.152C>A
        .variation("RYR1", "rs142474192", "G", "A")  // c.418G>A
    ;
    Path vcfFile = testWrapper.execute(null);

    testWrapper.testNotCalledByMatcher("CACNA1S");
    testWrapper.testCalledByMatcher("RYR1");

    List<String> ryr1ExpectedCalls = List.of("c.97A>G/[c.97A>G + c.152C>A + c.418G>A]");
    testWrapper.testSourceDiplotypes(DataSource.CPIC, "RYR1", ryr1ExpectedCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "RYR1", List.of("c.97A>G", "c.97A>G"));

    Document document = readHtmlReport(vcfFile);
    SortedMap<String, List<String>> expectedCallsMap = new TreeMap<>();
    expectedCallsMap.put("CACNA1S", NO_DATA);
    expectedCallsMap.put("RYR1", ryr1ExpectedCalls);
    htmlChecks(document, expectedCallsMap, "enflurane", RecPresence.YES, RecPresence.NO);
  }

  /**
   * Test a missing CACNA1S gene with 3 RYR1 het calls.
   */
  @Test
  void cacna1sMissing_ryr1Het3(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false, false, false);
    testWrapper.getVcfBuilder()
        .reference("RYR1")
        .variation("RYR1", "rs193922746", "A", "G") // c.97A>G
        .variation("RYR1", "rs193922749", "C", "A") // c.152C>A
        .variation("RYR1", "rs142474192", "G", "A")  // c.418G>A
    ;
    Path vcfFile = testWrapper.execute(null);

    testWrapper.testNotCalledByMatcher("CACNA1S");
    testWrapper.testCalledByMatcher("RYR1");

    List<String> ryr1ExpectedCalls = List.of("c.97A>G", "c.152C>A", "c.418G>A");
    testWrapper.testSourceDiplotypes(DataSource.CPIC, "RYR1", ryr1ExpectedCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "RYR1", List.of("c.97A>G", "c.152C>A"));

    Document document = readHtmlReport(vcfFile);
    SortedMap<String, List<String>> expectedCallsMap = new TreeMap<>();
    expectedCallsMap.put("CACNA1S", NO_DATA);
    expectedCallsMap.put("RYR1", ryr1ExpectedCalls);
    htmlChecks(document, expectedCallsMap, "enflurane", RecPresence.YES, RecPresence.NO);
  }

  /**
   * Test a missing CACNA1S gene with 3 RYR1 het calls.
   * Make sure malignant function variants are used to look up recommendations.
   */
  @Test
  void cacna1sMissing_ryr1Het4(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false, false, false);
    testWrapper.getVcfBuilder()
        .reference("RYR1")
        .variation("RYR1", "rs193922746", "A", "G") // c.97A>G - malignant
        .variation("RYR1", "rs193922749", "C", "A") // c.152C>A
        .variation("RYR1", "rs142474192", "G", "A") // c.418G>A
        .variation("RYR1", "rs146876145", "C", "T") // c.14918C>T - malignant
    ;
    Path vcfFile = testWrapper.execute(null);

    testWrapper.testNotCalledByMatcher("CACNA1S");
    testWrapper.testCalledByMatcher("RYR1");

    List<String> ryr1ExpectedCalls = List.of("c.97A>G", "c.152C>A", "c.418G>A", "c.14918C>T");
    testWrapper.testSourceDiplotypes(DataSource.CPIC, "RYR1", ryr1ExpectedCalls);
    // make sure malignant function variants are used to look up recommendations
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "RYR1", List.of("c.97A>G", "c.14918C>T"));

    Document document = readHtmlReport(vcfFile);
    SortedMap<String, List<String>> expectedCallsMap = new TreeMap<>();
    expectedCallsMap.put("CACNA1S", NO_DATA);
    expectedCallsMap.put("RYR1", ryr1ExpectedCalls);
    htmlChecks(document, expectedCallsMap, "enflurane", RecPresence.YES, RecPresence.NO);
  }

  /**
   * Test a missing CACNA1S gene with 3 RYR1 het calls, phased.
   * Make sure malignant function variants are used to look up recommendations.
   */
  @Test
  void cacna1sMissing_ryr1Het4Phased(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false, false, false);
    testWrapper.getVcfBuilder()
        .phased()
        .reference("RYR1")
        .variation("RYR1", "rs193922746", "A", "G") // c.97A>G - malignant
        .variation("RYR1", "rs193922749", "C", "A") // c.152C>A
        .variation("RYR1", "rs142474192", "G", "A") // c.418G>A
        .variation("RYR1", "rs146876145", "C", "T") // c.14918C>T - malignant
    ;
    Path vcfFile = testWrapper.execute(null);

    testWrapper.testNotCalledByMatcher("CACNA1S");
    testWrapper.testCalledByMatcher("RYR1");

    List<String> ryr1ExpectedCalls = List.of("Reference/[c.97A>G + c.152C>A + c.418G>A + c.14918C>T]");
    List<String> cpicStyleCalls = List.of("[c.97A>G + c.152C>A + c.418G>A + c.14918C>T] (heterozygous)");
    testWrapper.testSourceDiplotypes(DataSource.CPIC, "RYR1", cpicStyleCalls);
    // make sure malignant function variants are used to look up recommendations
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "RYR1", List.of("Reference", "c.97A>G"));

    Document document = readHtmlReport(vcfFile);
    SortedMap<String, List<String>> expectedCallsMap = new TreeMap<>();
    expectedCallsMap.put("CACNA1S", NO_DATA);
    expectedCallsMap.put("RYR1", ryr1ExpectedCalls);
    SortedMap<String, List<String>> cpicStyleCallsMap = new TreeMap<>();
    cpicStyleCallsMap.put("CACNA1S", NO_DATA);
    cpicStyleCallsMap.put("RYR1", cpicStyleCalls);
    htmlChecks(document, expectedCallsMap, cpicStyleCallsMap, "enflurane", RecPresence.YES, RecPresence.NO);
  }
}

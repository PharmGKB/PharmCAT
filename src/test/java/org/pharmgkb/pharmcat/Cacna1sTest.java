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
import org.pharmgkb.pharmcat.reporter.model.PrescribingGuidanceSource;

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
    TestUtils.setSaveTestOutput(true);
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
    Path vcfFile = testWrapper.execute();

    testWrapper.testCalledByMatcher("CACNA1S");
    testWrapper.testCalledByMatcher("RYR1");

    List<String> cacna1sExpectedCalls = List.of(TextConstants.REFERENCE + TextConstants.GENOTYPE_DELIMITER + TextConstants.REFERENCE);
    testWrapper.testSourceDiplotypes(DataSource.CPIC, "CACNA1S", cacna1sExpectedCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "CACNA1S", expectedCallsToRecommendedDiplotypes(cacna1sExpectedCalls));

    List<String> ryr1ExpectedCalls = List.of("c.97A>G", "c.152C>A");
    testWrapper.testSourceDiplotypes(DataSource.CPIC, "RYR1", ryr1ExpectedCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "RYR1", ryr1ExpectedCalls);

    Document document = readHtmlReport(vcfFile);
    htmlChecks(document,
        new ImmutableSortedMap.Builder<String, List<String>>(Ordering.natural())
            .put("CACNA1S", cacna1sExpectedCalls)
            .put("RYR1", ryr1ExpectedCalls)
            .build(),
        null,
        "enflurane", RecPresence.YES, RecPresence.NO);
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
    Path vcfFile = testWrapper.execute();

    testWrapper.testNotCalledByMatcher("RYR1");
    testWrapper.testCalledByMatcher("CACNA1S");

    List<String> cacna1sExpectedCalls = List.of(TextConstants.REFERENCE + TextConstants.GENOTYPE_DELIMITER + "c.520C>T");
    List<String> cacna1sCpicStyleCalls = List.of("c.520C>T (heterozygous)");

    testWrapper.testSourceDiplotypes(DataSource.CPIC, "CACNA1S", cacna1sCpicStyleCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "CACNA1S", expectedCallsToRecommendedDiplotypes(cacna1sExpectedCalls));
    testWrapper.testMatchedAnnotations("enflurane", PrescribingGuidanceSource.CPIC_GUIDELINE, 1);
    testWrapper.testNoMatchFromSource("desflurane", PrescribingGuidanceSource.DPWG_GUIDELINE);

    Document document = readHtmlReport(vcfFile);
    htmlChecks(document,
        new ImmutableSortedMap.Builder<String, List<String>>(Ordering.natural())
            .put("CACNA1S", cacna1sExpectedCalls)
            .put("RYR1", NO_DATA)
            .build(),
        new ImmutableSortedMap.Builder<String, List<String>>(Ordering.natural())
            .put("CACNA1S", cacna1sCpicStyleCalls)
            .build(),
         "enflurane", RecPresence.YES,
        new ImmutableSortedMap.Builder<String, String>(Ordering.natural())
            .put("CACNA1S", "Malignant Hyperthermia Susceptibility")
            .put("RYR1", "No Result")
            .build(),
        RecPresence.NO, null);
  }

  @Test
  void cacna1sHet2_ryr1Missing(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("CACNA1S", "rs1800559", "C", "T")
        .variation("CACNA1S", "rs772226819", "G", "A");
    Path vcfFile = testWrapper.execute();

    List<String> expectedCalls = List.of("c.520C>T/c.3257G>A");

    testWrapper.testCalledByMatcher("CACNA1S");
    testWrapper.testSourceDiplotypes(DataSource.CPIC, "CACNA1S", expectedCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "CACNA1S", expectedCallsToRecommendedDiplotypes(expectedCalls));
    testWrapper.testPrintCalls(DataSource.CPIC, "CACNA1S", expectedCalls);

    testWrapper.testMatchedAnnotations("desflurane", PrescribingGuidanceSource.CPIC_GUIDELINE, 1);
    testWrapper.testNoMatchFromSource("desflurane", PrescribingGuidanceSource.DPWG_GUIDELINE);

    Document document = readHtmlReport(vcfFile);
    htmlChecks(document,
        new ImmutableSortedMap.Builder<String, List<String>>(Ordering.natural())
            .put("CACNA1S", expectedCalls)
            .put("RYR1", NO_DATA)
            .build(),
        null, "desflurane", RecPresence.YES,
        new ImmutableSortedMap.Builder<String, String>(Ordering.natural())
            .put("CACNA1S", "Malignant Hyperthermia Susceptibility")
            .put("RYR1", "No Result")
            .build(),
        RecPresence.NO, null);
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
    Path vcfFile = testWrapper.execute();

    testWrapper.testCalledByMatcher("CACNA1S");
    testWrapper.testCalledByMatcher("RYR1");

    List<String> cacna1sExpectedCalls = List.of(TextConstants.REFERENCE + TextConstants.GENOTYPE_DELIMITER + "c.520C>T");
    List<String> cacna1sCpicStyleCalls = List.of("c.520C>T (heterozygous)");
    testWrapper.testSourceDiplotypes(DataSource.CPIC, "CACNA1S", cacna1sCpicStyleCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "CACNA1S", expectedCallsToRecommendedDiplotypes(cacna1sExpectedCalls));

    List<String> ryr1ExpectedCalls = List.of("Reference/c.7522C>G");
    List<String> ryr1CpicStyleCalls = List.of("c.7522C>G (heterozygous)");
    testWrapper.testSourceDiplotypes(DataSource.CPIC, "RYR1", ryr1CpicStyleCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "RYR1", expectedCallsToRecommendedDiplotypes(ryr1ExpectedCalls));

    Document document = readHtmlReport(vcfFile);
    htmlChecks(document,
        new ImmutableSortedMap.Builder<String, List<String>>(Ordering.natural())
            .put("CACNA1S", cacna1sExpectedCalls)
            .put("RYR1", ryr1ExpectedCalls)
            .build(),
        new ImmutableSortedMap.Builder<String, List<String>>(Ordering.natural())
            .put("CACNA1S", cacna1sCpicStyleCalls)
            .put("RYR1", ryr1CpicStyleCalls)
            .build(),
        "desflurane", RecPresence.YES,
        new ImmutableSortedMap.Builder<String, String>(Ordering.natural())
            .put("CACNA1S", "Malignant Hyperthermia Susceptibility")
            .put("RYR1", "Uncertain Susceptibility")
            .build(),
        RecPresence.NO, null);
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
    Path vcfFile = testWrapper.execute();

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
    SortedMap<String, List<String>> expectedCallsMap = new TreeMap<>();
    expectedCallsMap.put("CACNA1S", cacna1sExpectedCalls);
    expectedCallsMap.put("RYR1", ryr1ExpectedCalls);
    SortedMap<String, List<String>> cpicStyleCallsMap = new TreeMap<>();
    cpicStyleCallsMap.put("CACNA1S", cacna1sCpicStyleCalls);
    htmlChecks(document, expectedCallsMap, cpicStyleCallsMap, "enflurane", RecPresence.YES, RecPresence.NO);
  }
}

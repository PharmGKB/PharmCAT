package org.pharmgkb.pharmcat;

import java.nio.file.Path;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Stream;
import com.google.common.base.Preconditions;
import org.checkerframework.checker.nullness.qual.Nullable;
import org.jsoup.nodes.Document;
import org.jsoup.select.Elements;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.TestInfo;
import org.pharmgkb.pharmcat.haplotype.model.CombinationMatch;
import org.pharmgkb.pharmcat.reporter.TextConstants;
import org.pharmgkb.pharmcat.reporter.model.DataSource;
import org.pharmgkb.pharmcat.reporter.model.result.AnnotationReport;
import org.pharmgkb.pharmcat.reporter.model.result.DrugReport;
import org.pharmgkb.pharmcat.reporter.model.result.GeneReport;
import org.pharmgkb.pharmcat.util.HaplotypeNameComparator;

import static org.junit.jupiter.api.Assertions.*;
import static org.pharmgkb.pharmcat.PipelineTest.*;


/**
 * This JUnit tests DPYD.
 *
 * @author Mark Woon
 */
class DpydTest {


  @BeforeAll
  static void prepare() {
    //TestUtils.setSaveTestOutput(true);
  }

  @AfterEach
  void deleteDirectory(TestInfo testInfo) {
    TestUtils.deleteTestOutputDirectory(testInfo);
  }


  private @Nullable List<String> callsToComponents(List<String> expectedCalls) {
    if (expectedCalls == null || expectedCalls.size() != 1) {
      return null;
    }
    if (expectedCalls.get(0).equals("Unknown/Unknown") ||
        !expectedCalls.get(0).contains(TextConstants.GENOTYPE_DELIMITER)) {
      return null;
    }
    return expectedCalls.stream()
        .flatMap(s -> Arrays.stream(s.split(TextConstants.GENOTYPE_DELIMITER)))
        .flatMap(s -> {
          if (s.startsWith("[")) {
            s = s.substring(1, s.length() - 1);
            return CombinationMatch.COMBINATION_NAME_SPLITTER.splitToList(s).stream();
          }
          return Stream.of(s);
        })
        .distinct()
        .sorted(new HaplotypeNameComparator())
        .toList();
  }


  private void dpydHasReports(PipelineWrapper testWrapper, RecPresence hasDpwgReport) {
    dpydHasReports(testWrapper, RecPresence.YES, hasDpwgReport);
  }

  private void dpydHasReports(PipelineWrapper testWrapper, RecPresence hasCpicReport, RecPresence hasDpwgReport) {
    GeneReport cpicDpydGeneReport = testWrapper.getContext().getGeneReport(DataSource.CPIC, "DPYD");
    assertNotNull(cpicDpydGeneReport);
    assertEquals(1, cpicDpydGeneReport.getRecommendationDiplotypes().size());

    GeneReport dpwgDpydGeneReport = testWrapper.getContext().getGeneReport(DataSource.DPWG, "DPYD");
    assertNotNull(dpwgDpydGeneReport);
    assertEquals(1, dpwgDpydGeneReport.getRecommendationDiplotypes().size());

    int numAnnotations = 0;

    if (hasCpicReport == RecPresence.YES) {
      testWrapper.testAnyMatchFromSource("fluorouracil", DataSource.CPIC);
      testWrapper.testAnyMatchFromSource("capecitabine", DataSource.CPIC);
      numAnnotations += 1;

    } else {
      testWrapper.testNoMatchFromSource("fluorouracil", DataSource.CPIC);
      testWrapper.testNoMatchFromSource("capecitabine", DataSource.CPIC);
    }

    if (hasDpwgReport == RecPresence.YES) {
      testWrapper.testAnyMatchFromSource("fluorouracil", DataSource.DPWG);
      testWrapper.testAnyMatchFromSource("capecitabine", DataSource.DPWG);
      numAnnotations += 1;

    } else {
      testWrapper.testNoMatchFromSource("fluorouracil", DataSource.DPWG);
      testWrapper.testNoMatchFromSource("capecitabine", DataSource.DPWG);
    }

    if (numAnnotations > 0) {
      testWrapper.testMatchedAnnotations("fluorouracil", numAnnotations);
      testWrapper.testMatchedAnnotations("capecitabine", numAnnotations);
    }
  }



  @Test
  void testDpydPhased(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .phased()
        .variation("DPYD", "rs3918290", "C", "T")
        .variation("DPYD", "rs1801159", "C", "T")
    ;
    Path vcfFile = testWrapper.execute(null);

    String gene = "DPYD";
    List<String> expectedCalls = List.of("c.1627A>G (*5)/c.1905+1G>A (*2A)");
    RecPresence hasDpwgAnnotations = RecPresence.YES;

    testWrapper.testCalledByMatcher(gene);
    testWrapper.testSourceDiplotypes(DataSource.CPIC, gene, expectedCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, gene, expectedCallsToRecommendedDiplotypes(expectedCalls));
    testWrapper.testPrintCalls(DataSource.CPIC, gene, expectedCalls);

    dpydHasReports(testWrapper, hasDpwgAnnotations);

    Document document = readHtmlReport(vcfFile);
    dpydHtmlChecks(document, expectedCalls, false, hasDpwgAnnotations);
  }


  @Test
  void testDpydUnphased(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("DPYD", "rs3918290", "C", "T")
        .variation("DPYD", "rs1801159", "C", "T");
    Path vcfFile = testWrapper.execute(null);

    List<String> expectedCalls = List.of("c.1627A>G (*5)", "c.1905+1G>A (*2A)");
    RecPresence hasDpwgAnnotations = RecPresence.YES;

    testWrapper.testCalledByMatcher("DPYD");
    testWrapper.testSourceDiplotypes(DataSource.CPIC, "DPYD", expectedCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "DPYD", expectedCalls);
    testWrapper.testPrintCalls(DataSource.CPIC, "DPYD", expectedCalls);

    dpydHasReports(testWrapper, hasDpwgAnnotations);

    Document document = readHtmlReport(vcfFile);
    dpydHtmlChecks(document, expectedCalls, false, hasDpwgAnnotations);
  }


  @Test
  void testDpydUndocumentedVariation(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .allowUnknownAllele()
        .variation("DPYD", "rs3918290", "C", "G")
    ;
    Path vcfFile = testWrapper.execute(null);

    List<String> expectedCalls = List.of("Reference/Reference");
    RecPresence hasDpwgAnnotations = RecPresence.YES;

    testWrapper.testCalledByMatcher("DPYD");
    testWrapper.testSourceDiplotypes(DataSource.CPIC, "DPYD", expectedCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "DPYD", expectedCallsToRecommendedDiplotypes(expectedCalls));
    testWrapper.testPrintCalls(DataSource.CPIC, "DPYD", expectedCalls);

    dpydHasReports(testWrapper, hasDpwgAnnotations);

    Document document = readHtmlReport(vcfFile);
    dpydHtmlChecks(document, expectedCalls, false, hasDpwgAnnotations);
  }


  private void dpydHtmlChecks(Document document, @Nullable List<String> expectedCalls, boolean hasMissingPositions,
      RecPresence hasDpwgAnnotation) {
    dpydHtmlChecks(document, expectedCalls, hasMissingPositions, RecPresence.YES, hasDpwgAnnotation);
  }

  /**
   * Checks for expected HTML output for DPYD.  Only drug checked is capecitabine.
   */
  private void dpydHtmlChecks(Document document, @Nullable List<String> expectedCalls, boolean hasMissingPositions,
      RecPresence hasCpicAnnotation, RecPresence hasDpwgAnnotation) {

    boolean noCall = expectedCalls != null && expectedCalls.size() == 1 &&
        expectedCalls.get(0).equals("Unknown/Unknown");
    List<String> expectedComponents = callsToComponents(expectedCalls);

    if (expectedComponents != null) {
      if (expectedCalls != null) {
        Elements gsLeastFunction = document.select(".gs-DPYD .gs-dip_leastFunction");
        assertEquals(expectedCalls.size(), gsLeastFunction.size());
        assertEquals(expectedCalls,
            gsLeastFunction.stream()
                .map(e -> e.child(0).text())
                .toList());
      }
      Elements gsComponents = document.select(".gs-DPYD .gs-dip_component");
      assertEquals(expectedComponents.size(), gsComponents.size());
      List<String> components = gsComponents.stream()
          .map(e -> e.child(0).text())
          .toList();
      assertEquals(expectedComponents, components);

    } else {
      Elements gsDips = document.select(".gs-DPYD .gs-dip");
      if (noCall) {
        assertEquals(0, gsDips.size());
      } else {
        Preconditions.checkNotNull(expectedCalls);
        assertEquals(expectedCalls.size(), gsDips.size());
        assertEquals(expectedCalls,
            gsDips.stream()
                .map(e -> e.child(0).text())
                .toList());
      }
    }

    Elements capecitabineSection = document.getElementsByClass("capecitabine");
    if (noCall) {
      assertEquals(0, capecitabineSection.size());
    } else {
      assertEquals(1, capecitabineSection.size());
      // should have DPYD warning
      Elements capecitabineMsgs = capecitabineSection.get(0).getElementsByClass("alert-info");
      assertEquals(hasMissingPositions ? 2 : 1, capecitabineMsgs.size());
      assertTrue(capecitabineMsgs.get(0).text().contains("lowest activity"));

      if (expectedCalls != null) {
        List<String> expectedRxCalls = expectedCalls.stream()
            .map(c -> "DPYD:" + c)
            .toList();

        Elements cpicCapecitabineDips = capecitabineSection.select(".cpic-capecitabine .rx-dip");
        if (hasCpicAnnotation == RecPresence.YES) {
          assertEquals(expectedRxCalls,
              cpicCapecitabineDips.stream()
                  .map(e -> cleanupRxDip(e, List.of("DPYD")))
                  .toList());

        } else {
          assertEquals(0, cpicCapecitabineDips.size());
          Elements unmatchedDips = capecitabineSection.select(".cpic-capecitabine .rx-unmatched-dip");
          assertEquals(expectedRxCalls, unmatchedDips.stream()
              .map(e -> cleanupRxDip(e, List.of("DPYD")))
              .toList());
        }

        Elements dpwgCapecitabineDips = capecitabineSection.select(".dpwg-capecitabine .rx-dip");
        if (hasDpwgAnnotation == RecPresence.YES) {
          assertEquals(expectedRxCalls,
              dpwgCapecitabineDips.stream()
                  .map(e -> cleanupRxDip(e, List.of("DPYD")))
                  .toList());
        } else {
          assertEquals(0, dpwgCapecitabineDips.size());
          Elements unmatchedDips = capecitabineSection.select(".dpwg-capecitabine .rx-unmatched-dip");
          assertEquals(expectedRxCalls, unmatchedDips.stream()
              .map(e -> cleanupRxDip(e, List.of("DPYD")))
              .toList());
        }
      }
    }

    Elements dpydSection = document.select(".gene.dpyd");
    assertEquals(1, dpydSection.size());
    assertEquals(0, dpydSection.get(0).getElementsByClass("no-data").size());
    Elements gsResult = dpydSection.select(".genotype-result");
    assertEquals(1, gsResult.size());
    if (noCall) {
      assertEquals(TextConstants.UNCALLED, gsResult.get(0).text());
    }
  }


  @Test
  void testDpydUnphasedMultiple1(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("DPYD", "rs183385770", "C", "T")  // 0 activity value
        .variation("DPYD", "rs186169810", "A", "C") // 0.5 activity value
        .variation("DPYD", "rs112766203", "G", "A") // c.2279C>T - 0.5 activity value
        .variation("DPYD", "rs144395748", "G", "C"); // 1 activity value
    Path vcfFile = testWrapper.execute(null);

    List<String> expectedCalls = List.of(
        "c.1024G>A",
        "c.1314T>G",
        "c.1358C>G",
        "c.2279C>T"
    );
    RecPresence hasDpwgAnnotations = RecPresence.NO;

    testWrapper.testCalledByMatcher("DPYD");
    testWrapper.testSourceDiplotypes(DataSource.CPIC, "DPYD", expectedCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "DPYD", List.of("c.1024G>A", "c.1314T>G"));
    testWrapper.testPrintCalls(DataSource.CPIC, "DPYD", expectedCalls);
    testWrapper.testLookupByActivity(DataSource.CPIC, "DPYD", "0.5");

    dpydHasReports(testWrapper, hasDpwgAnnotations);

    Document document = readHtmlReport(vcfFile);
    dpydHtmlChecks(document, expectedCalls, false, hasDpwgAnnotations);
  }

  @Test
  void testDpydUnphasedMultiple2(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("DPYD", "rs183385770", "C", "T")  // 0 activity value
        .variation("DPYD", "rs186169810", "A", "C") // 0.5 activity value
        .variation("DPYD", "rs67376798", "T", "A") // c.2846A>T
        .variation("DPYD", "rs144395748", "G", "C"); // 1 activity value

    Path vcfFile = testWrapper.execute(null);

    List<String> expectedCalls = List.of(
        "c.1024G>A",
        "c.1314T>G",
        "c.1358C>G",
        "c.2846A>T"
    );
    RecPresence hasDpwgAnnotations = RecPresence.NO;

    testWrapper.testCalledByMatcher("DPYD");
    testWrapper.testSourceDiplotypes(DataSource.CPIC, "DPYD", expectedCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "DPYD", List.of("c.1024G>A", "c.2846A>T"));
    testWrapper.testPrintCalls(DataSource.CPIC, "DPYD", expectedCalls);
    testWrapper.testLookupByActivity(DataSource.CPIC, "DPYD", "0.5");

    dpydHasReports(testWrapper, hasDpwgAnnotations);

    Document document = readHtmlReport(vcfFile);
    dpydHtmlChecks(document, expectedCalls, false, hasDpwgAnnotations);
  }

  @Test
  void testDpydC2846het(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .phased()
        .variation("DPYD", "rs67376798", "T", "A");
    Path vcfFile = testWrapper.execute(null);

    List<String> expectedCalls = List.of("Reference/c.2846A>T");
    RecPresence hasDpwgAnnotations = RecPresence.YES;

    testWrapper.testCalledByMatcher("DPYD");
    testWrapper.testSourceDiplotypes(DataSource.CPIC, "DPYD", expectedCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "DPYD", expectedCallsToRecommendedDiplotypes(expectedCalls));
    testWrapper.testPrintCalls(DataSource.CPIC, "DPYD", expectedCalls);

    dpydHasReports(testWrapper, hasDpwgAnnotations);

    Document document = readHtmlReport(vcfFile);
    dpydHtmlChecks(document, expectedCalls, false, hasDpwgAnnotations);
  }

  /**
   * Test to make sure AS 1.0 and AS 1.5 have the same recommendations for capecitabine. This is due to an update given
   * to the guideline after publication. See the
   * <a href="https://cpicpgx.org/guidelines/guideline-for-fluoropyrimidines-and-dpyd/">November 2018 update</a> for
   * details.
   */
  @Test
  void testDpydDifferenceOnScore(TestInfo testInfo) throws Exception {
    PipelineWrapper highScoreWrapper = new PipelineWrapper(testInfo, false);
    highScoreWrapper.getVcfBuilder()
        .phased()
        .variation("DPYD", "rs67376798", "T", "A");
    highScoreWrapper.execute(null);

    highScoreWrapper.testCalledByMatcher("DPYD");
    highScoreWrapper.testPrintCpicCalls("DPYD", "Reference/c.2846A>T");
    highScoreWrapper.testRecommendedDiplotypes("DPYD", "Reference", "c.2846A>T");
    GeneReport highScoreDpydReport = highScoreWrapper.getContext().getGeneReport(DataSource.CPIC, "DPYD");
    assertNotNull(highScoreDpydReport);
    assertTrue(highScoreDpydReport.getRecommendationDiplotypes().stream().allMatch((d) -> d.getActivityScore().equals("1.5")));

    highScoreWrapper.testAnyMatchFromSource("capecitabine", DataSource.CPIC);
    DrugReport highScoreDrug = highScoreWrapper.getContext().getDrugReport(DataSource.CPIC, "capecitabine");
    assertNotNull(highScoreDrug);
    List<String> highRecs = highScoreDrug.getGuidelines().stream()
        .flatMap(g -> g.getAnnotations().stream())
        .map(AnnotationReport::getDrugRecommendation)
        .toList();
    assertEquals(1, highRecs.size());

    PipelineWrapper lowScoreWrapper = new PipelineWrapper(testInfo, false);
    lowScoreWrapper.getVcfBuilder()
        .phased()
        .variation("DPYD", "rs67376798", "A", "A");
    lowScoreWrapper.execute(null);
    lowScoreWrapper.testCalledByMatcher("DPYD");
    lowScoreWrapper.testPrintCpicCalls("DPYD", "c.2846A>T/c.2846A>T");
    lowScoreWrapper.testRecommendedDiplotypes("DPYD", "c.2846A>T", "c.2846A>T");
    GeneReport lowScoreDpydReport = lowScoreWrapper.getContext().getGeneReport(DataSource.CPIC, "DPYD");
    assertNotNull(lowScoreDpydReport);
    assertTrue(lowScoreDpydReport.getRecommendationDiplotypes().stream().allMatch((d) -> d.getActivityScore().equals("1.0")));

    lowScoreWrapper.testAnyMatchFromSource("capecitabine", DataSource.CPIC);
    DrugReport lowScoreDrug = lowScoreWrapper.getContext().getDrugReport(DataSource.CPIC, "capecitabine");
    assertNotNull(lowScoreDrug);
    List<String> lowRecs = lowScoreDrug.getGuidelines().stream()
        .flatMap(g -> g.getAnnotations().stream())
        .map(AnnotationReport::getDrugRecommendation)
        .toList();
    assertEquals(1, lowRecs.size());

    // this is the point of this test:
    assertTrue(lowRecs.get(0).startsWith(highRecs.get(0)));
  }

  /**
   * This test puts 2 alleles on each strand of a phased DPYD and then asserts that the least-function allele is used
   * for lookup on each of the strands.
   */
  @Test
  void testDpydPhasedMultiTrans(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .phased()
        .variation("DPYD", "rs67376798", "A", "T") // Strand 1 decreased - c.2846A>T
        .variation("DPYD", "rs72547601", "C", "T") // Strand 1 no function - c.2933A>G
        .variation("DPYD", "rs60139309", "T", "C") // Strand 2 normal function - c.2582A>G
        .variation("DPYD", "rs139834141", "C", "T") // Strand 2 normal function - c.498G>A
    ;
    Path vcfFile = testWrapper.execute(null);

    List<String> expectedCalls = List.of("[c.498G>A + c.2582A>G]/[c.2846A>T + c.2933A>G]");
    RecPresence hasDpwgAnnotations = RecPresence.NO;

    testWrapper.testCalledByMatcher("DPYD");
    testWrapper.testSourceDiplotypes(DataSource.CPIC, "DPYD", expectedCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "DPYD", List.of("c.498G>A", "c.2933A>G"));
    testWrapper.testPrintCalls(DataSource.CPIC, "DPYD", expectedCalls);

    dpydHasReports(testWrapper, hasDpwgAnnotations);

    Document document = readHtmlReport(vcfFile);
    dpydHtmlChecks(document, expectedCalls, false, hasDpwgAnnotations);
  }

  /**
   * This test is the same as the previous test but DPYD is unphased instead of phased. This means the individual found
   * alleles should be reported and then the two least-function alleles should be used for recommendation lookup.
   */
  @Test
  void testDpydUnphasedMultiTrans(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("DPYD", "rs67376798", "A", "T") // decreased - c.2846A>T
        .variation("DPYD", "rs72547601", "C", "T") // no function - c.2933A>G
        .variation("DPYD", "rs60139309", "T", "C") // normal function - c.2582A>G
        .variation("DPYD", "rs139834141", "C", "T") // normal function - c.498G>A
    ;
    Path vcfFile = testWrapper.execute(null);

    List<String> expectedCalls = List.of("c.498G>A", "c.2582A>G", "c.2846A>T", "c.2933A>G");
    RecPresence hasDpwgAnnotations = RecPresence.NO;

    testWrapper.testCalledByMatcher("DPYD");
    testWrapper.testSourceDiplotypes(DataSource.CPIC, "DPYD", expectedCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "DPYD", List.of("c.2933A>G", "c.2846A>T"));
    testWrapper.testPrintCalls(DataSource.CPIC, "DPYD", expectedCalls);

    dpydHasReports(testWrapper, hasDpwgAnnotations);

    Document document = readHtmlReport(vcfFile);
    dpydHtmlChecks(document, expectedCalls, false, hasDpwgAnnotations);
  }

  @Test
  void testDpydS12het(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .phased()
        .variation("DPYD", "rs78060119", "C", "A");
    Path vcfFile = testWrapper.execute(null);

    List<String> expectedCalls = List.of("Reference/c.1156G>T (*12)");
    RecPresence hasDpwgAnnotations = RecPresence.NO;

    testWrapper.testCalledByMatcher("DPYD");
    testWrapper.testSourceDiplotypes(DataSource.CPIC, "DPYD", expectedCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "DPYD", expectedCallsToRecommendedDiplotypes(expectedCalls));
    testWrapper.testPrintCalls(DataSource.CPIC, "DPYD", expectedCalls);

    dpydHasReports(testWrapper, hasDpwgAnnotations);

    Document document = readHtmlReport(vcfFile);
    dpydHtmlChecks(document, expectedCalls, false, hasDpwgAnnotations);
  }

  @Test
  void testDpydHomNoFunctionEffectivelyPhased(TestInfo testInfo) throws Exception {
    // effectively phased
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("DPYD", "rs72549310", "A", "A")   // c.61C>T, hom variant (No function)
        .variation("DPYD", "rs150385342", "C", "T"); // c.313G>A het variant (Normal function)
    Path vcfFile = testWrapper.execute(null);

    List<String> expectedCalls = List.of("c.61C>T/[c.61C>T + c.313G>A]");
    RecPresence hasDpwgAnnotations = RecPresence.NO;

    testWrapper.testCalledByMatcher("DPYD");
    testWrapper.testSourceDiplotypes(DataSource.CPIC, "DPYD", expectedCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "DPYD", List.of("c.61C>T", "c.61C>T"));
    testWrapper.testPrintCalls(DataSource.CPIC, "DPYD", expectedCalls);

    dpydHasReports(testWrapper, hasDpwgAnnotations);

    Document document = readHtmlReport(vcfFile);
    dpydHtmlChecks(document, expectedCalls, false, hasDpwgAnnotations);
  }

  @Test
  void testDpydHomNoFunctionPhased(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .phased()
        .variation("DPYD", "rs72549310", "A", "A")   // c.61C>T, hom variant (No function)
        .variation("DPYD", "rs150385342", "C", "T"); // c.313G>A het variant (Normal function)
    Path vcfFile = testWrapper.execute(null);

    List<String> expectedCalls = List.of("c.61C>T/[c.61C>T + c.313G>A]");
    RecPresence hasDpwgAnnotations = RecPresence.NO;

    testWrapper.testCalledByMatcher("DPYD");
    testWrapper.testSourceDiplotypes(DataSource.CPIC, "DPYD", expectedCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "DPYD", List.of("c.61C>T", "c.61C>T"));
    testWrapper.testPrintCalls(DataSource.CPIC, "DPYD", expectedCalls);

    dpydHasReports(testWrapper, hasDpwgAnnotations);

    Document document = readHtmlReport(vcfFile);
    dpydHtmlChecks(document, expectedCalls, false, hasDpwgAnnotations);
  }

  @Test
  void testDpydHomNoFunctionUnphased(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("DPYD", "rs72547601", "C", "C") // c.2933A>G - no function
        .variation("DPYD", "rs67376798", "A", "T") // c.2846A>T - decreased
        .variation("DPYD", "rs60139309", "T", "C") // c.2582A>G - normal
    ;
    Path vcfFile = testWrapper.execute(null);

    List<String> expectedCalls = List.of("c.2582A>G", "c.2846A>T", "c.2933A>G");
    RecPresence hasDpwgAnnotations = RecPresence.NO;

    testWrapper.testCalledByMatcher("DPYD");
    testWrapper.testSourceDiplotypes(DataSource.CPIC, "DPYD", expectedCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "DPYD", List.of("c.2933A>G", "c.2933A>G"));
    testWrapper.testPrintCalls(DataSource.CPIC, "DPYD", expectedCalls);

    dpydHasReports(testWrapper, hasDpwgAnnotations);

    Document document = readHtmlReport(vcfFile);
    dpydHtmlChecks(document, expectedCalls, false, hasDpwgAnnotations);
  }

  @Test
  void testDpydHapB3_het_wobble_het(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("DPYD", "rs56038477", "C", "T") // g.97573863C>T
        .variation("DPYD", "rs75017182", "G", "C") // g.97579893G>C
    ;
    Path vcfFile = testWrapper.execute(null);

    List<String> expectedCalls = List.of("Reference", "c.1129-5923C>G, c.1236G>A (HapB3)");
    RecPresence hasDpwgAnnotations = RecPresence.YES;

    testWrapper.testCalledByMatcher("DPYD");
    testWrapper.testSourceDiplotypes(DataSource.CPIC, "DPYD", expectedCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "DPYD", expectedCalls);
    testWrapper.testPrintCalls(DataSource.CPIC, "DPYD", expectedCalls);

    dpydHasReports(testWrapper, hasDpwgAnnotations);

    Document document = readHtmlReport(vcfFile);
    dpydHtmlChecks(document, expectedCalls, false, hasDpwgAnnotations);
  }

  @Test
  void testDpydHapB3_het_wobble_C_hom(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("DPYD", "rs56038477", "C", "C") // g.97573863C>T
        .variation("DPYD", "rs75017182", "G", "C") // g.97579893G>C
    ;
    Path vcfFile = testWrapper.execute(null);

    List<String> expectedCalls = List.of("Reference/c.1129-5923C>G, c.1236G>A (HapB3)");
    RecPresence hasDpwgAnnotations = RecPresence.YES;

    testWrapper.testCalledByMatcher("DPYD");
    testWrapper.testSourceDiplotypes(DataSource.CPIC, "DPYD", expectedCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "DPYD", expectedCallsToRecommendedDiplotypes(expectedCalls));
    testWrapper.testPrintCalls(DataSource.CPIC, "DPYD", expectedCalls);

    dpydHasReports(testWrapper, hasDpwgAnnotations);

    Document document = readHtmlReport(vcfFile);
    dpydHtmlChecks(document, expectedCalls, false, hasDpwgAnnotations);
  }


  @Test
  void testDpydHapB3_hom_wobble_T_hom(TestInfo testInfo) throws Exception {
    // effectively phased, homozygous alternative
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("DPYD", "rs56038477", "T", "T") // g.97573863C>T
        .variation("DPYD", "rs75017182", "C", "C") // g.97579893G>C
    ;
    Path vcfFile = testWrapper.execute(null);

    List<String> expectedCalls = List.of("c.1129-5923C>G, c.1236G>A (HapB3)/c.1129-5923C>G, c.1236G>A (HapB3)");
    RecPresence hasDpwgAnnotations = RecPresence.YES;

    testWrapper.testCalledByMatcher("DPYD");
    testWrapper.testSourceDiplotypes(DataSource.CPIC, "DPYD", expectedCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "DPYD", expectedCallsToRecommendedDiplotypes(expectedCalls));
    testWrapper.testPrintCalls(DataSource.CPIC, "DPYD", expectedCalls);

    dpydHasReports(testWrapper, hasDpwgAnnotations);

    Document document = readHtmlReport(vcfFile);
    dpydHtmlChecks(document, expectedCalls, false, hasDpwgAnnotations);
  }

  @Test
  void testDpydHapB3_hom_wobble_C_hom(TestInfo testInfo) throws Exception {
    // effectively phased, homozygous alternative
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("DPYD", "rs56038477", "C", "C") // g.97573863C>T
        .variation("DPYD", "rs75017182", "C", "C") // g.97579893G>C
    ;
    Path vcfFile = testWrapper.execute(null);

    List<String> expectedCalls = List.of("c.1129-5923C>G, c.1236G>A (HapB3)/c.1129-5923C>G, c.1236G>A (HapB3)");
    RecPresence hasDpwgAnnotations = RecPresence.YES;

    testWrapper.testCalledByMatcher("DPYD");
    testWrapper.testSourceDiplotypes(DataSource.CPIC, "DPYD", expectedCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "DPYD", expectedCallsToRecommendedDiplotypes(expectedCalls));
    testWrapper.testPrintCalls(DataSource.CPIC, "DPYD", expectedCalls);

    dpydHasReports(testWrapper, hasDpwgAnnotations);

    Document document = readHtmlReport(vcfFile);
    dpydHtmlChecks(document, expectedCalls, false, hasDpwgAnnotations);
  }

  @Test
  void testDpydHapB3_hom_wobble_het(TestInfo testInfo) throws Exception {
    // effectively phased, homozygous alternative
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("DPYD", "rs56038477", "C", "T") // g.97573863C>T
        .variation("DPYD", "rs75017182", "C", "C") // g.97579893G>C
    ;
    Path vcfFile = testWrapper.execute(null);

    List<String> expectedCalls = List.of("c.1129-5923C>G, c.1236G>A (HapB3)/c.1129-5923C>G, c.1236G>A (HapB3)");
    RecPresence hasDpwgAnnotations = RecPresence.YES;

    testWrapper.testCalledByMatcher("DPYD");
    testWrapper.testSourceDiplotypes(DataSource.CPIC, "DPYD", expectedCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "DPYD", expectedCallsToRecommendedDiplotypes(expectedCalls));
    testWrapper.testPrintCalls(DataSource.CPIC, "DPYD", expectedCalls);

    dpydHasReports(testWrapper, hasDpwgAnnotations);

    Document document = readHtmlReport(vcfFile);
    dpydHtmlChecks(document, expectedCalls, false, hasDpwgAnnotations);
  }

  @Test
  void testDpydHapB3_only_wobble_C_hom_gets_reference(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("DPYD", "rs56038477", "C", "C") // g.97573863C>T
    //.variation("DPYD", "rs75017182", "G", "C") // g.97579893G>C
    ;
    Path vcfFile = testWrapper.execute(null);

    List<String> expectedCalls = List.of("Reference/Reference");
    RecPresence hasDpwgAnnotations = RecPresence.YES;

    testWrapper.testCalledByMatcher("DPYD");
    testWrapper.testSourceDiplotypes(DataSource.CPIC, "DPYD", expectedCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "DPYD", expectedCallsToRecommendedDiplotypes(expectedCalls));
    testWrapper.testPrintCalls(DataSource.CPIC, "DPYD", expectedCalls);

    dpydHasReports(testWrapper, hasDpwgAnnotations);

    Document document = readHtmlReport(vcfFile);
    dpydHtmlChecks(document, expectedCalls, false, hasDpwgAnnotations);
  }

  @Test
  void testDpydHapB3_noCall_only_wobble_T_hom(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("DPYD", "rs56038477", "T", "T") // g.97573863C>T
    //.variation("DPYD", "rs75017182", "G", "C") // g.97579893G>C
    ;
    Path vcfFile = testWrapper.execute(null);

    List<String> expectedCalls = List.of("Unknown/Unknown");
    RecPresence hasDpwgAnnotations = RecPresence.NO;

    testWrapper.testNotCalledByMatcher("DPYD");
    testWrapper.testSourceDiplotypes(DataSource.CPIC, "DPYD", expectedCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "DPYD", expectedCalls);
    testWrapper.testPrintCalls(DataSource.CPIC, "DPYD", List.of(TextConstants.UNCALLED));

    dpydHasReports(testWrapper, RecPresence.NO, hasDpwgAnnotations);

    Document document = readHtmlReport(vcfFile);
    dpydHtmlChecks(document, expectedCalls, false, RecPresence.NO, hasDpwgAnnotations);
  }

  @Test
  void testDpydHapB3_noCall_only_wobble_het(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("DPYD", "rs56038477", "C", "T") // g.97573863C>T
    //.variation("DPYD", "rs75017182", "G", "C") // g.97579893G>C
    ;
    Path vcfFile = testWrapper.execute(null);

    List<String> expectedCalls = List.of("Unknown/Unknown");
    RecPresence hasDpwgAnnotations = RecPresence.NO;

    testWrapper.testNotCalledByMatcher("DPYD");
    testWrapper.testSourceDiplotypes(DataSource.CPIC, "DPYD", expectedCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "DPYD", expectedCalls);
    testWrapper.testPrintCalls(DataSource.CPIC, "DPYD", List.of(TextConstants.UNCALLED));

    dpydHasReports(testWrapper, RecPresence.NO, hasDpwgAnnotations);

    Document document = readHtmlReport(vcfFile);
    dpydHtmlChecks(document, expectedCalls, false, RecPresence.NO, hasDpwgAnnotations);
  }


  @Test
  void testDpydHapB3_wobble_het_nonRef_hom(TestInfo testInfo) throws Exception {
    TestUtils.setSaveTestOutput(true);
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        // c.61C>T named allele
        .variation("DPYD", "rs72549310", "A", "A")   // c.61C>T, hom variant (No function)
        // partial HapB3
        .variation("DPYD", "rs56038477", "C", "T") // g.97573863C>T
    //.variation("DPYD", "rs75017182", "G", "C") // g.97579893G>C
    ;
    Path vcfFile = testWrapper.execute(null);

    List<String> expectedCalls = List.of("c.61C>T/c.61C>T");
    RecPresence hasDpwgAnnotations = RecPresence.NO;

    testWrapper.testCalledByMatcher("DPYD");
    testWrapper.testSourceDiplotypes(DataSource.CPIC, "DPYD", expectedCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "DPYD", expectedCallsToRecommendedDiplotypes(expectedCalls));
    testWrapper.testPrintCalls(DataSource.CPIC, "DPYD", expectedCalls);

    dpydHasReports(testWrapper, hasDpwgAnnotations);

    Document document = readHtmlReport(vcfFile);
    dpydHtmlChecks(document, expectedCalls, false, hasDpwgAnnotations);
  }

  @Test
  void testDpydHapB3_wobble_het_nonRef(TestInfo testInfo) throws Exception {
    TestUtils.setSaveTestOutput(true);
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("DPYD", "rs72549310", "G", "A")   // c.61C>T, hom variant (No function)
        .variation("DPYD", "rs56038477", "C", "T")   // g.97573863C>T
    //.variation("DPYD", "rs75017182", "G", "C") // g.97579893G>C
    ;
    Path vcfFile = testWrapper.execute(null);

    List<String> expectedCalls = List.of("c.61C>T");
    RecPresence hasCpicAnnotations = RecPresence.YES_NO_MATCH;
    RecPresence hasDpwgAnnotations = RecPresence.YES_NO_MATCH;

    testWrapper.testCalledByMatcher("DPYD");
    testWrapper.testSourceDiplotypes(DataSource.CPIC, "DPYD", expectedCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "DPYD", expectedCalls);
    testWrapper.testPrintCalls(DataSource.CPIC, "DPYD", expectedCalls);

    dpydHasReports(testWrapper, hasCpicAnnotations, hasDpwgAnnotations);

    Document document = readHtmlReport(vcfFile);
    dpydHtmlChecks(document, expectedCalls, false, hasCpicAnnotations, hasDpwgAnnotations);
  }

  @Test
  void testDpydHapB3_rs75017182_missing(TestInfo testInfo) throws Exception {
    // effectively phased
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, true, false, false);
    testWrapper.getVcfBuilder()
        .variation("DPYD", "rs56038477", "C", "T") // g.97573863C>T
        .missing("DPYD", "rs75017182")
    ;
    Path vcfFile = testWrapper.execute(null);

    List<String> expectedCalls = List.of("Reference", "c.1129-5923C>G, c.1236G>A (HapB3)");
    RecPresence hasDpwgAnnotations = RecPresence.YES;

    testWrapper.testCalledByMatcher("DPYD");
    testWrapper.testSourceDiplotypes(DataSource.CPIC, "DPYD", expectedCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "DPYD", expectedCalls);
    testWrapper.testPrintCalls(DataSource.CPIC, "DPYD", expectedCalls);

    dpydHasReports(testWrapper, hasDpwgAnnotations);

    Document document = readHtmlReport(vcfFile);
    dpydHtmlChecks(document, expectedCalls, true, hasDpwgAnnotations);
  }

  @Test
  void testDpydHapB3_rs56038477_missing(TestInfo testInfo) throws Exception {
    // effectively phased
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("DPYD", "rs75017182", "G", "C") // g.97579893G>C
        .missing("DPYD", "rs56038477")
    ;
    Path vcfFile = testWrapper.execute(null);

    List<String> expectedCalls = List.of("Reference/c.1129-5923C>G, c.1236G>A (HapB3)");
    RecPresence hasDpwgAnnotations = RecPresence.YES;

    testWrapper.testCalledByMatcher("DPYD");
    testWrapper.testSourceDiplotypes(DataSource.CPIC, "DPYD", expectedCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "DPYD", expectedCallsToRecommendedDiplotypes(expectedCalls));
    testWrapper.testPrintCalls(DataSource.CPIC, "DPYD", expectedCalls);

    dpydHasReports(testWrapper, hasDpwgAnnotations);

    Document document = readHtmlReport(vcfFile);
    dpydHtmlChecks(document, expectedCalls, true, hasDpwgAnnotations);
  }
}

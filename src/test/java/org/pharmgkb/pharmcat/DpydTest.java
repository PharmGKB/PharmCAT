package org.pharmgkb.pharmcat;

import java.io.BufferedReader;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.List;
import java.util.Objects;
import java.util.stream.Stream;
import com.google.common.base.Preconditions;
import org.apache.commons.lang3.StringUtils;
import org.checkerframework.checker.nullness.qual.Nullable;
import org.jsoup.nodes.Document;
import org.jsoup.nodes.Element;
import org.jsoup.select.Elements;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.TestInfo;
import org.pharmgkb.common.util.PathUtils;
import org.pharmgkb.pharmcat.definition.model.VariantLocus;
import org.pharmgkb.pharmcat.haplotype.DpydHapB3Matcher;
import org.pharmgkb.pharmcat.haplotype.model.CombinationMatch;
import org.pharmgkb.pharmcat.reporter.MessageHelper;
import org.pharmgkb.pharmcat.reporter.TextConstants;
import org.pharmgkb.pharmcat.reporter.handlebars.ReportHelpers;
import org.pharmgkb.pharmcat.reporter.model.DataSource;
import org.pharmgkb.pharmcat.reporter.model.PrescribingGuidanceSource;
import org.pharmgkb.pharmcat.reporter.model.VariantReport;
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
    ReportHelpers.setDebugMode(true);
    //TestUtils.setSaveTestOutput(true);
  }

  @AfterEach
  void deleteDirectory(TestInfo testInfo) {
    TestUtils.deleteTestOutputDirectory(testInfo);
  }


  public static @Nullable List<String> callsToComponents(List<String> expectedCalls) {
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


  private static void dpydHasReports(PipelineWrapper testWrapper, RecPresence hasDpwgReport) {
    dpydHasReports(testWrapper, RecPresence.YES, hasDpwgReport);
  }

  static void dpydHasReports(PipelineWrapper testWrapper, RecPresence hasCpicReport, RecPresence hasDpwgReport) {
    GeneReport cpicDpydGeneReport = testWrapper.getContext().getGeneReport(DataSource.CPIC, "DPYD");
    assertNotNull(cpicDpydGeneReport);
    assertEquals(1, cpicDpydGeneReport.getRecommendationDiplotypes().size());

    GeneReport dpwgDpydGeneReport = testWrapper.getContext().getGeneReport(DataSource.DPWG, "DPYD");
    assertNotNull(dpwgDpydGeneReport);
    assertEquals(1, dpwgDpydGeneReport.getRecommendationDiplotypes().size());

    if (hasCpicReport == RecPresence.YES) {
      testWrapper.testAnyMatchFromSource("fluorouracil", PrescribingGuidanceSource.CPIC_GUIDELINE);
      testWrapper.testAnyMatchFromSource("capecitabine", PrescribingGuidanceSource.CPIC_GUIDELINE);
      testWrapper.testMatchedAnnotations("fluorouracil", PrescribingGuidanceSource.CPIC_GUIDELINE, 1);
      testWrapper.testMatchedAnnotations("capecitabine", PrescribingGuidanceSource.CPIC_GUIDELINE, 1);

    } else {
      testWrapper.testNoMatchFromSource("fluorouracil", PrescribingGuidanceSource.CPIC_GUIDELINE);
      testWrapper.testNoMatchFromSource("capecitabine", PrescribingGuidanceSource.CPIC_GUIDELINE);
    }

    if (hasDpwgReport == RecPresence.YES) {
      testWrapper.testAnyMatchFromSource("fluorouracil", PrescribingGuidanceSource.DPWG_GUIDELINE);
      testWrapper.testAnyMatchFromSource("capecitabine", PrescribingGuidanceSource.DPWG_GUIDELINE);
      testWrapper.testMatchedAnnotations("fluorouracil", PrescribingGuidanceSource.DPWG_GUIDELINE, 1);
      testWrapper.testMatchedAnnotations("capecitabine", PrescribingGuidanceSource.DPWG_GUIDELINE, 1);

    } else {
      testWrapper.testNoMatchFromSource("fluorouracil", PrescribingGuidanceSource.DPWG_GUIDELINE);
      testWrapper.testNoMatchFromSource("capecitabine", PrescribingGuidanceSource.DPWG_GUIDELINE);
    }
  }


  private void doStandardChecks(PipelineWrapper testWrapper, Path vcfFile, List<String> expectedCalls,
      boolean hasMissingPositions) throws Exception {
    doStandardChecks(testWrapper, vcfFile, expectedCalls, null, null, hasMissingPositions,
        RecPresence.YES);
  }

  private void doStandardChecks(PipelineWrapper testWrapper, Path vcfFile, List<String> expectedCalls,
      @Nullable List<String> cpicStyleCalls, @Nullable List<String> recommendedDips, boolean hasMissingPositions,
      RecPresence hasDpwgAnnotations) throws Exception {

    if (cpicStyleCalls == null) {
      cpicStyleCalls = expectedCalls;
    }
    if (recommendedDips == null) {
      if (expectedCalls.size() == 1) {
        recommendedDips = expectedCallsToRecommendedDiplotypes(expectedCalls);
      } else {
        recommendedDips = expectedCalls;
      }
    }
    testWrapper.testCalledByMatcher("DPYD");
    testWrapper.testSourceDiplotypes(DataSource.CPIC, "DPYD", expectedCalls, cpicStyleCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "DPYD", recommendedDips);
    testWrapper.testPrintCalls(DataSource.CPIC, "DPYD", cpicStyleCalls);

    dpydHasReports(testWrapper, hasDpwgAnnotations);

    Document document = readHtmlReport(vcfFile);
    dpydHtmlChecks(document, expectedCalls, cpicStyleCalls, hasMissingPositions, RecPresence.YES, hasDpwgAnnotations);
  }


  @Test
  void testDpydPhased(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .phased()
        .variation("DPYD", "rs3918290", "C", "T")
        .variation("DPYD", "rs1801159", "C", "T")
    ;
    Path vcfFile = testWrapper.execute();

    List<String> expectedCalls = List.of("c.1627A>G (*5)/c.1905+1G>A (*2A)");

    doStandardChecks(testWrapper, vcfFile, expectedCalls, false);

    GeneReport dpwgReport = testWrapper.getContext().getGeneReport(DataSource.DPWG, "DPYD");
    assertNotNull(dpwgReport);
    assertTrue(dpwgReport.getRecommendationDiplotypes().stream().flatMap((d) -> d.getLookupKeys().stream()).noneMatch(TextConstants::isUnspecified), "DPWG missing lookup key for DPYD");

    // DPWG does not include the 1627 allele in function definition so use Reference for lookup
    String gene = "DPYD";
    testWrapper.testRecommendedDiplotypes(DataSource.DPWG, gene, expectedCallsToRecommendedDiplotypes(List.of("Reference/c.1905+1G>A (*2A)")));
    // all other diplotype usage can use the alleles as called
    testWrapper.testSourceDiplotypes(DataSource.DPWG, gene, expectedCalls);
    testWrapper.testPrintCalls(DataSource.DPWG, gene, expectedCalls);
  }


  @Test
  void testDpydUnphased(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("DPYD", "rs3918290", "C", "T")
        .variation("DPYD", "rs1801159", "C", "T");
    Path vcfFile = testWrapper.execute();

    List<String> expectedCalls = List.of("c.1627A>G (*5)", "c.1905+1G>A (*2A)");

    doStandardChecks(testWrapper, vcfFile, expectedCalls, false);
  }


  @Test
  void testDpydUndocumentedVariation(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .allowUnknownAllele()
        .variation("DPYD", "rs3918290", "C", "G")
    ;
    Path vcfFile = testWrapper.execute();

    List<String> expectedCalls = List.of("Reference/Reference");

    doStandardChecks(testWrapper, vcfFile, expectedCalls, false);
  }


  /**
   * Checks for expected HTML output for DPYD.
   * The only drug checked is capecitabine.
   */
  static void dpydHtmlChecks(Document document, @Nullable List<String> expectedCalls,
      @Nullable List<String> cpicStyleCalls, boolean hasMissingPositions,
      RecPresence hasCpicAnnotation, RecPresence hasDpwgAnnotation) {

    boolean noCall = expectedCalls != null && expectedCalls.size() == 1 &&
        expectedCalls.get(0).equals("Unknown/Unknown");
    List<String> expectedComponents = callsToComponents(expectedCalls);

    if (expectedComponents != null) {
      Elements gsLowestFunction = document.select(".gs-DPYD .gs-dip_lowestFunction");
      assertEquals(expectedCalls.size(), gsLowestFunction.size());
      assertEquals(cpicStyleCalls == null ? expectedCalls : cpicStyleCalls,
          gsLowestFunction.stream()
              .map(e -> e.child(0).text())
              .toList());

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
      Elements missingVariantsWarning = capecitabineSection.get(0).select(".alert-info.missing-variants");
      if (hasMissingPositions) {
        assertEquals(1, missingVariantsWarning.size(), "Expected missing variants warning, but did not find it");
      } else {
        assertEquals(0, missingVariantsWarning.size(), "Not expecting missing variants warning, but found it");
      }

      List<String> expectedRxCalls = cpicStyleCalls == null ?
          expectedCalls.stream()
              .map(c -> "DPYD:" + c)
              .toList() :
          cpicStyleCalls.stream()
              .map(c -> "DPYD:" + c)
              .toList();

      Elements cpicCapecitabineDips = capecitabineSection.select(".cpic-guideline-capecitabine .rx-dip");
      if (hasCpicAnnotation == RecPresence.YES) {
        assertEquals(expectedRxCalls,
            cpicCapecitabineDips.stream()
                .map(e -> cleanupRxDip(e, List.of("DPYD")))
                .toList());

      } else {
        assertEquals(0, cpicCapecitabineDips.size());
        Elements unmatchedDips = capecitabineSection.select(".cpic-guideline-capecitabine .rx-unmatched-dip");
        assertEquals(expectedRxCalls, unmatchedDips.stream()
            .map(e -> cleanupRxDip(e, List.of("DPYD")))
            .toList());
      }

      Elements dpwgCapecitabineDips = capecitabineSection.select(".dpwg-guideline-capecitabine .rx-dip");
      if (hasDpwgAnnotation == RecPresence.YES) {
        assertEquals(expectedRxCalls,
            dpwgCapecitabineDips.stream()
                .map(e -> cleanupRxDip(e, List.of("DPYD")))
                .toList());
      } else {
        assertEquals(0, dpwgCapecitabineDips.size());
        Elements unmatchedDips = capecitabineSection.select(".dpwg-guideline-capecitabine .rx-unmatched-dip");
        assertEquals(expectedRxCalls, unmatchedDips.stream()
            .map(e -> cleanupRxDip(e, List.of("DPYD")))
            .toList());
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
    Path vcfFile = testWrapper.execute();

    List<String> expectedCalls = List.of(
        "c.1024G>A",
        "c.1314T>G",
        "c.1358C>G",
        "c.2279C>T"
    );

    doStandardChecks(testWrapper, vcfFile, expectedCalls, null, List.of("c.1024G>A", "c.1314T>G"), false, RecPresence.NO);
    testWrapper.testLookupByActivity(DataSource.CPIC, "DPYD", "0.5");
  }

  @Test
  void testDpydUnphasedMultiple2(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("DPYD", "rs183385770", "C", "T")  // 0 activity value
        .variation("DPYD", "rs186169810", "A", "C") // 0.5 activity value
        .variation("DPYD", "rs67376798", "T", "A") // c.2846A>T
        .variation("DPYD", "rs144395748", "G", "C"); // 1 activity value

    Path vcfFile = testWrapper.execute();

    List<String> expectedCalls = List.of(
        "c.1024G>A",
        "c.1314T>G",
        "c.1358C>G",
        "c.2846A>T"
    );

    doStandardChecks(testWrapper, vcfFile, expectedCalls, null, List.of("c.1024G>A", "c.2846A>T"),
        false, RecPresence.NO);

    testWrapper.testLookupByActivity(DataSource.CPIC, "DPYD", "0.5");

    // in this instance, tegafur should have only a DPWG annotation, but it has no matching guidance
    DrugReport tegafur = testWrapper.getContext().getDrugReport(PrescribingGuidanceSource.DPWG_GUIDELINE, "tegafur");
    assertNotNull(tegafur);
    assertFalse(tegafur.getGuidelines().stream().anyMatch(g -> g.getSource() == PrescribingGuidanceSource.CPIC_GUIDELINE));
    assertTrue(tegafur.getGuidelines().stream().anyMatch(g -> g.getSource() == PrescribingGuidanceSource.DPWG_GUIDELINE));
    assertFalse(tegafur.getGuidelines().stream().anyMatch(g -> g.getSource() == PrescribingGuidanceSource.FDA_LABEL));
    assertFalse(tegafur.getGuidelines().stream().anyMatch(g -> g.getSource() == PrescribingGuidanceSource.FDA_ASSOC));
    assertFalse(tegafur.isMatched());
  }

  @Test
  void testDpydC2846het(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .phased()
        .variation("DPYD", "rs67376798", "T", "A");
    Path vcfFile = testWrapper.execute();

    List<String> expectedCalls = List.of("Reference/c.2846A>T");
    List<String> cpicStyleCalls = List.of("c.2846A>T (heterozygous)");

    doStandardChecks(testWrapper, vcfFile, expectedCalls, cpicStyleCalls, null, false, RecPresence.YES);
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
    highScoreWrapper.execute();

    highScoreWrapper.testCalledByMatcher("DPYD");
    highScoreWrapper.testPrintCpicCalls("DPYD", "c.2846A>T (heterozygous)");
    highScoreWrapper.testRecommendedDiplotypes("DPYD", "Reference", "c.2846A>T");
    GeneReport highScoreDpydReport = highScoreWrapper.getContext().getGeneReport(DataSource.CPIC, "DPYD");
    assertNotNull(highScoreDpydReport);
    assertTrue(highScoreDpydReport.getRecommendationDiplotypes().stream().allMatch((d) -> d.getActivityScore().equals("1.5")));

    highScoreWrapper.testAnyMatchFromSource("capecitabine", PrescribingGuidanceSource.CPIC_GUIDELINE);
    DrugReport highScoreDrug = highScoreWrapper.getContext().getDrugReport(PrescribingGuidanceSource.CPIC_GUIDELINE, "capecitabine");
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
    lowScoreWrapper.execute();
    lowScoreWrapper.testCalledByMatcher("DPYD");
    lowScoreWrapper.testPrintCpicCalls("DPYD", "c.2846A>T/c.2846A>T");
    lowScoreWrapper.testRecommendedDiplotypes("DPYD", "c.2846A>T", "c.2846A>T");
    GeneReport lowScoreDpydReport = lowScoreWrapper.getContext().getGeneReport(DataSource.CPIC, "DPYD");
    assertNotNull(lowScoreDpydReport);
    assertTrue(lowScoreDpydReport.getRecommendationDiplotypes().stream().allMatch((d) -> d.getActivityScore().equals("1.0")));

    lowScoreWrapper.testAnyMatchFromSource("capecitabine", PrescribingGuidanceSource.CPIC_GUIDELINE);
    DrugReport lowScoreDrug = lowScoreWrapper.getContext().getDrugReport(PrescribingGuidanceSource.CPIC_GUIDELINE, "capecitabine");
    assertNotNull(lowScoreDrug);
    List<String> lowRecs = lowScoreDrug.getGuidelines().stream()
        .flatMap(g -> g.getAnnotations().stream())
        .map(AnnotationReport::getDrugRecommendation)
        .toList();
    assertEquals(1, lowRecs.size());

    // this is the point of this test:
    String lowRec = lowRecs.get(0);
    String highRec = highRecs.get(0);
    // we're doing this split because there is "other consideration" text in this field that IS different between them
    assertEquals(lowRec.split("\\(if available\\)")[0], highRec.split("\\(if available\\)")[0]);
  }

  /**
   * This test puts 2 alleles on each strand of a phased DPYD and then asserts that the
   * lowest-function allele is used for lookup on each of the strands.
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
    Path vcfFile = testWrapper.execute();

    List<String> expectedCalls = List.of("[c.498G>A + c.2582A>G]/[c.2846A>T + c.2933A>G]");

    doStandardChecks(testWrapper, vcfFile, expectedCalls, null, List.of("c.498G>A", "c.2933A>G"), false, RecPresence.NO);
  }

  /**
   * This test is the same as the previous test, but DPYD is unphased instead of phased.
   * This means the individual found alleles should be reported and then the two lowest-function
   * alleles should be used for recommendation lookup.
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
    Path vcfFile = testWrapper.execute();

    List<String> expectedCalls = List.of("c.498G>A", "c.2582A>G", "c.2846A>T", "c.2933A>G");

    doStandardChecks(testWrapper, vcfFile, expectedCalls, null, List.of("c.2933A>G", "c.2846A>T"), false, RecPresence.NO);
  }

  @Test
  void testDpydUnphasedMultiTrans_infer_homo(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("DPYD", "rs67376798", "A", "T") // decreased - c.2846A>T
        .variation("DPYD", "rs72547601", "C", "C") // no function - c.2933A>G
        .variation("DPYD", "rs60139309", "T", "C") // normal function - c.2582A>G
        .variation("DPYD", "rs139834141", "C", "T") // normal function - c.498G>A
    ;
    Path vcfFile = testWrapper.execute();

    List<String> expectedCalls = List.of("c.498G>A", "c.2582A>G", "c.2846A>T", "c.2933A>G (homozygous)");

    doStandardChecks(testWrapper, vcfFile, expectedCalls, null, List.of("c.2933A>G", "c.2933A>G"), false, RecPresence.NO);
  }

  @Test
  void testDpydS12het(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .phased()
        .variation("DPYD", "rs78060119", "C", "A");
    Path vcfFile = testWrapper.execute();

    List<String> expectedCalls = List.of("Reference/c.1156G>T (*12)");
    List<String> cpicStyleCalls = List.of("c.1156G>T (*12) (heterozygous)");

    doStandardChecks(testWrapper, vcfFile, expectedCalls, cpicStyleCalls, null, false, RecPresence.NO);
  }

  @Test
  void testDpydHomNoFunctionEffectivelyPhased(TestInfo testInfo) throws Exception {
    // effectively phased
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("DPYD", "rs72549310", "A", "A")   // c.61C>T, hom variant (No function)
        .variation("DPYD", "rs150385342", "C", "T"); // c.313G>A het variant (Normal function)
    Path vcfFile = testWrapper.execute();

    List<String> expectedCalls = List.of("c.61C>T/[c.61C>T + c.313G>A]");

    doStandardChecks(testWrapper, vcfFile, expectedCalls, null, List.of("c.61C>T", "c.61C>T"), false, RecPresence.NO);
  }

  @Test
  void testDpydHomNoFunctionPhased(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .phased()
        .variation("DPYD", "rs72549310", "A", "A")   // c.61C>T, hom variant (No function)
        .variation("DPYD", "rs150385342", "C", "T"); // c.313G>A het variant (Normal function)
    Path vcfFile = testWrapper.execute();

    List<String> expectedCalls = List.of("c.61C>T/[c.61C>T + c.313G>A]");

    doStandardChecks(testWrapper, vcfFile, expectedCalls, null, List.of("c.61C>T", "c.61C>T"), false, RecPresence.NO);
  }

  @Test
  void testDpydHomNoFunctionUnphased(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("DPYD", "rs72547601", "C", "C") // c.2933A>G - no function
        .variation("DPYD", "rs67376798", "A", "T") // c.2846A>T - decreased
        .variation("DPYD", "rs60139309", "T", "C") // c.2582A>G - normal
    ;
    Path vcfFile = testWrapper.execute();

    List<String> expectedCalls = List.of("c.2582A>G", "c.2846A>T", "c.2933A>G (homozygous)");

    doStandardChecks(testWrapper, vcfFile, expectedCalls, null, List.of("c.2933A>G", "c.2933A>G"),
        false, RecPresence.NO);
  }


  @Test
  void test155(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .phased()
        // chr1	97515865	1|0  C->T (*4) 1.0
        .variation("DPYD", "rs1801158", "T", "C")
        // chr1	97573863	0|1 C->T
        .variation("DPYD", "rs56038477", "C", "T")
        // chr1	97579893	0|1 G->C
        .variation("DPYD", "rs75017182", "G", "C")
        // chr1	97699535	1|0 T->C (c.496A>G) 1.0
        .variation("DPYD", "rs2297595", "C", "T")
        // chr1	97883329	1|1 A->G (*9A) 1.0
        .variation("DPYD", "rs1801265", "G", "G")
    ;

    Path vcfFile = testWrapper.execute();

    List<String> expectedCalls = List.of(
        "[c.85T>C (*9A) + c.1129-5923C>G, c.1236G>A (HapB3)]/[c.85T>C (*9A) + c.496A>G + c.1601G>A (*4)]"
    );
    List<String> recommendedDiplotypes = List.of("c.85T>C (*9A)", DpydHapB3Matcher.HAPB3_ALLELE);

    doStandardChecks(testWrapper, vcfFile, expectedCalls, null, recommendedDiplotypes, false, RecPresence.YES);
  }

  @Test
  void test156(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .phased()
        // chr1	97305364	1|0  C->T
        .variation("DPYD", "rs1801160", "T", "C")
        // chr1	97515839	1|0 T->C
        .variation("DPYD", "rs1801159", "C", "T")

        // chr1	97573863	1|0 C->T
        .variation("DPYD", "rs56038477", "T", "C")
        // chr1	97579893	1|0 G->C
        .variation("DPYD", "rs75017182", "C", "G")
    ;

    Path vcfFile = testWrapper.execute();

    List<String> expectedCalls = List.of(
        "Reference/[c.1129-5923C>G, c.1236G>A (HapB3) + c.1627A>G (*5) + c.2194G>A (*6)]"
    );
    List<String> cpicStyleCalls = List.of(
        "[c.1129-5923C>G, c.1236G>A (HapB3) + c.1627A>G (*5) + c.2194G>A (*6)] (heterozygous)"
    );

    doStandardChecks(testWrapper, vcfFile, expectedCalls, cpicStyleCalls,
        List.of("Reference", DpydHapB3Matcher.HAPB3_ALLELE), false, RecPresence.YES);
  }


  /**
   * Check that DPYD matching isn't affected when the "find-combinations" mode is enabled.
   */
  @Test
  void testFindCombinations(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, true, false, false);
    testWrapper.getVcfBuilder()
        .phased()
        .allowUnknownAllele()
        .variation("DPYD", "rs114096998", "G", "A")
    ;

    Path vcfFile = testWrapper.execute();

    List<String> expectedCalls = List.of(
        "Reference/Reference"
    );

    doStandardChecks(testWrapper, vcfFile, expectedCalls, false);
  }


  @Test
  void hapB3AndIntronicC(TestInfo testInfo) throws Exception {

    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, true, false, false);
    testWrapper.getVcfBuilder()
        .phased()
        // hapB3 exon C>T
        .variation("DPYD", "rs56038477", "C", "T")
        // hapB3 intron G>C
        .variation("DPYD", "rs75017182", "C", "C")
    ;

    Path vcfFile = testWrapper.execute();

    List<String> expectedCalls = List.of(
        "c.1129-5923C>G/c.1129-5923C>G, c.1236G>A (HapB3)"
    );
    doStandardChecks(testWrapper, vcfFile, expectedCalls, false);
  }

  @Test
  void hapB3AndIntronicC_missingFirstAllele(TestInfo testInfo) throws Exception {

    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, true, false, false)
        .saveIntermediateFiles();
    testWrapper.getVcfBuilder()
        .phased()
        // hapB3 exon C>T
        .variation("DPYD", "rs56038477", ".", "T")
        // hapB3 intron G>C
        .variation("DPYD", "rs75017182", ".", "C")
    ;

    Path vcfFile = testWrapper.execute();

    List<String> expectedCalls = List.of(
        "Reference/Reference"
    );

    doStandardChecks(testWrapper, vcfFile, expectedCalls, true);
  }

  @Test
  void hapB3AndIntronicC_missingSecondAllele(TestInfo testInfo) throws Exception {

    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, true, false, false)
        .saveIntermediateFiles();
    testWrapper.getVcfBuilder()
        .phased()
        // hapB3 exon C>T
        .variation("DPYD", "rs56038477", "T", ".")
        // hapB3 intron G>C
        .variation("DPYD", "rs75017182", "C", ".")
    ;

    Path vcfFile = testWrapper.execute();

    List<String> expectedCalls = List.of(
        "Reference/c.1129-5923C>G, c.1236G>A (HapB3)"
    );
    List<String> cpicStyleCalls = List.of(
        "c.1129-5923C>G, c.1236G>A (HapB3) (heterozygous)"
    );
    doStandardChecks(testWrapper, vcfFile, expectedCalls, cpicStyleCalls, null, true, RecPresence.YES);
  }

  @Test
  void hapB3AndIntronicC_missingAlternateStrands1(TestInfo testInfo) throws Exception {

    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, true, false, false)
        .saveIntermediateFiles();
    testWrapper.getVcfBuilder()
        .phased()
        // hapB3 exon C>T
        .variation("DPYD", "rs56038477", "T", ".")
        // hapB3 intron G>C
        .variation("DPYD", "rs75017182", ".", "C")
    ;

    Path vcfFile = testWrapper.execute();

    List<String> expectedCalls = List.of(
        "Reference/c.1129-5923C>G"
    );
    List<String> cpicStyleCalls = List.of(
        "c.1129-5923C>G (heterozygous)"
    );

    doStandardChecks(testWrapper, vcfFile, expectedCalls, cpicStyleCalls, null, true, RecPresence.YES);
  }

  @Test
  void hapB3AndIntronicC_missingAlternateStrands2(TestInfo testInfo) throws Exception {

    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, true, false, false)
        .saveIntermediateFiles();
    testWrapper.getVcfBuilder()
        .phased()
        // hapB3 exon C>T
        .variation("DPYD", "rs56038477", ".", "T")
        // hapB3 intron G>C
        .variation("DPYD", "rs75017182", "C", ".")
    ;

    Path vcfFile = testWrapper.execute();

    List<String> expectedCalls = List.of(
        "Reference/c.1129-5923C>G"
    );
    List<String> cpicStyleCalls = List.of(
        "c.1129-5923C>G (heterozygous)"
    );

    doStandardChecks(testWrapper, vcfFile, expectedCalls, cpicStyleCalls, null, true, RecPresence.YES);
  }


  @Test
  void hapB3AndIntronicD(TestInfo testInfo) throws Exception {

    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, true, false, false);
    testWrapper.getVcfBuilder()
        .phased()
        // hapB3 exon C>T
        .variation("DPYD", "rs56038477", "T", "C")
        // hapB3 intron G>C
        .variation("DPYD", "rs75017182", "C", "C")
    ;

    Path vcfFile = testWrapper.execute();

    List<String> expectedCalls = List.of(
        "c.1129-5923C>G/c.1129-5923C>G, c.1236G>A (HapB3)"
    );
    List<String> recommendedDiplotypes = List.of("c.1129-5923C>G, c.1236G>A (HapB3)", "c.1129-5923C>G");
    doStandardChecks(testWrapper, vcfFile, expectedCalls, null, recommendedDiplotypes, false, RecPresence.YES);
  }


  @Test
  void hapB3(TestInfo testInfo) throws Exception {

    // this file based on https://docs.google.com/spreadsheets/d/1M3XDbrCmgz7RCZFMvWgn7wYjiFufvi9JPtDhhqh_Bc8/
    Path tsvFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/DPYD/hapb3.tsv");
    try (BufferedReader reader = Files.newBufferedReader(tsvFile)) {
      // header row
      String line = reader.readLine();
      assertNotNull(line);
      // name row
      line = reader.readLine();
      assertNotNull(line);
      // rsid row
      line = reader.readLine();
      assertNotNull(line);
      String[] data = line.split("\t");
      assertEquals("rs56038477", data[0]);
      assertEquals("rs75017182", data[1]);
      assertEquals("rs72549310", data[2]);
      assertEquals("rs67376798", data[3]);
      assertEquals("rs1801265", data[4]);
      assertEquals("rs186169810", data[5]);
      assertEquals("rs115232898", data[6]);
      int row = 3;
      while ((line = reader.readLine()) != null) {
        row += 1;
        if (row < 99) {
          //continue;
        }
        String name = String.format("row_%03d", row);
        data = line.split("\t");

        System.out.println("\n" + name);
        assertTrue(data.length > 14, name + " only has " + data.length + " columns, expecting at least 15");

        runDpydTest(testInfo, name, data, true);
        runDpydTest(testInfo, name, data, false);
      }
    }
  }

  private void runDpydTest(TestInfo testInfo, String name, String[] data, boolean isPhased) throws Exception {
    if (isPhased) {
      name += " - phased";
      System.out.println("\tphased");
    } else {
      name += " - unphased";
      System.out.println("\tunphased");
    }
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, name, false, false, false);
    TestVcfBuilder vcfBuilder = testWrapper.getVcfBuilder();
    if (isPhased) {
      vcfBuilder.phased();
    }
    setVariation(testWrapper, vcfBuilder, "rs56038477", data[0]);
    setVariation(testWrapper, vcfBuilder, "rs75017182", data[1]);
    setVariation(testWrapper, vcfBuilder, "rs72549310", data[2]);
    setVariation(testWrapper, vcfBuilder, "rs67376798", data[3]);
    setVariation(testWrapper, vcfBuilder, "rs1801265", data[4]);
    setVariation(testWrapper, vcfBuilder, "rs186169810", data[5], "G");
    setVariation(testWrapper, vcfBuilder, "rs115232898", data[6]);
    assertFalse(data[7].isEmpty());

    Path vcfFile = testWrapper.execute();

    // 7 - phased call
    // 8 - phased recommendation
    // 9 - unphased call
    // 10 - unphased recommendation

    List<String> expectedCalls = parseDiplotypes(data[7]);
    if (!isPhased && !data[9].isEmpty()) {
      expectedCalls = parseDiplotypes(data[9]);
    }
    List<String> recommendedDiplotypes;
    if (isPhased) {
      if (!data[8].isEmpty()) {
        recommendedDiplotypes = parseDiplotypes(data[8]);
      } else {
        recommendedDiplotypes = expectedCalls;
      }
    } else {
      // not phased
      if (!data[10].isEmpty()) {
        recommendedDiplotypes = parseDiplotypes(data[10]);
      } else if (!data[8].isEmpty()) {
        recommendedDiplotypes = parseDiplotypes(data[8]);
      } else {
        recommendedDiplotypes = expectedCalls;
      }
    }
    if (recommendedDiplotypes.size() == 1) {
      recommendedDiplotypes = expectedCallsToRecommendedDiplotypes(recommendedDiplotypes);
    }

    List<String> cpicStyleCalls = expectedCalls.stream()
        .map(d -> {
          if (d.startsWith("Reference/") && !d.endsWith("/Reference")) {
            return d.substring(10) + " (heterozygous)";
          }
          return null;
        })
        .filter(Objects::nonNull)
        .toList();
    if (cpicStyleCalls.isEmpty()) {
      cpicStyleCalls = null;
    }

    testWrapper.testCalledByMatcher("DPYD");
    testWrapper.testSourceDiplotypes(DataSource.CPIC, "DPYD", expectedCalls, cpicStyleCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "DPYD", recommendedDiplotypes);
    testWrapper.testPrintCalls(DataSource.CPIC, "DPYD", cpicStyleCalls == null ? expectedCalls : cpicStyleCalls);

    GeneReport cpicGeneReport = testWrapper.getContext().getGeneReport(DataSource.CPIC, "DPYD");
    assertNotNull(cpicGeneReport);
    int numMissing = 0;
    if ("missing".equals(data[0])) {
      numMissing += 1;
    }
    if ("missing".equals(data[1])) {
      numMissing += 1;
    }
    if ("missing".equals(data[3])) {
      numMissing += 1;
    }
    if ("missing".equals(data[4])) {
      numMissing += 1;
    }
    if ("missing".equals(data[5])) {
      numMissing += 1;
    }
    if ("missing".equals(data[6])) {
      numMissing += 1;
    }
    long numMissingSampleAlleles = cpicGeneReport.getVariantReports().stream()
        .filter(VariantReport::isMissing)
        .count();
    assertEquals(numMissing, numMissingSampleAlleles);

    Document document = readHtmlReport(vcfFile);
    Element dpydSection = document.select(".gene.dpyd").get(0);

    // check warnings
    // 11 - pgx warning
    // 12 - phased warning
    // 13 - unphased warning

    // pgx warning
    if ("yes".equalsIgnoreCase(data[11])) {
      assertTrue(cpicGeneReport.getVariantReports().stream()
          .anyMatch(vr -> !vr.getWarnings().isEmpty()), "Expecting PGx warnings, found none");
    } else {
      assertTrue(cpicGeneReport.getVariantReports().stream()
          .allMatch(vr -> vr.getWarnings().isEmpty()), "Should not have PGx warnings, but found them");
    }

    String warning = isPhased ? StringUtils.stripToNull(data[12]) : StringUtils.stripToNull(data[13]);
    if (warning != null) {
      assertTrue(cpicGeneReport.getMessages().stream()
          .anyMatch(ma -> ma.getName().equals(warning)), "Missing '" + warning + "' warning");
      assertEquals(1, dpydSection.getElementsByClass(warning).size());
    } else {
      List<String> warnings = cpicGeneReport.getMessages().stream()
          .map(ma -> {
            if (ma.getName().equals(MessageHelper.MSG_DPYD_HAPB3_EXONIC_ONLY) ||
                ma.getName().equals(MessageHelper.MSG_DPYD_HAPB3_INTRONIC_MISMATCH_EXONIC)) {
              return ma.getName();
            }
            return null;
          })
          .filter(Objects::nonNull)
          .toList();

      assertTrue(warnings.isEmpty(), "Expected no warnings, found " + String.join(", ", warnings));
      assertEquals(0, dpydSection.getElementsByClass(MessageHelper.MSG_DPYD_HAPB3_EXONIC_ONLY).size());
      assertEquals(0, dpydSection.getElementsByClass(MessageHelper.MSG_DPYD_HAPB3_INTRONIC_MISMATCH_EXONIC).size());
    }

    RecPresence hasCpicAnnotations = RecPresence.fromString(data[14]);
    RecPresence hasDpwgAnnotations = RecPresence.fromString(data[15]);
    dpydHasReports(testWrapper, hasCpicAnnotations, hasDpwgAnnotations);
    dpydHtmlChecks(document, expectedCalls, cpicStyleCalls, numMissingSampleAlleles > 0, hasCpicAnnotations, hasDpwgAnnotations);
  }


  private static List<String> parseDiplotypes(String data) {
    return Arrays.stream(data.split(","))
        .map(a -> StringUtils.stripToEmpty(a)
            .replaceAll("(\\S)\\+(\\S)", "$1 + $2")
            .replaceAll("ref", "Reference")
            .replaceAll("HapB3Intron", DpydHapB3Matcher.HAPB3_INTRONIC_ALLELE)
            .replaceAll("HapB3", DpydHapB3Matcher.HAPB3_ALLELE))
        .toList();
  }


  private void setVariation(PipelineWrapper testWrapper, TestVcfBuilder vcfBuilder, String rsid, String callValue) {
    setVariation(testWrapper, vcfBuilder, rsid, callValue, null);
  }

  private void setVariation(PipelineWrapper testWrapper, TestVcfBuilder vcfBuilder, String rsid, String callValue,
      @Nullable String nonPgxAllele) {

    if (callValue.equalsIgnoreCase("missing")) {
      vcfBuilder.missing("DPYD", rsid);
      return;
    }
    String[] calls;
    if (callValue.contains("/")) {
      calls = callValue.split("/");
    } else if (callValue.contains("|")) {
      calls = callValue.split("\\|");
    } else {
      throw new IllegalArgumentException("Invalid call: " + callValue);
    }
    assertEquals(2, calls.length);

    VariantLocus vl = Arrays.stream(testWrapper.getEnv().getDefinitionReader().getPositions("DPYD"))
        .filter(v -> v.getRsid().equals(rsid))
        .findAny()
        .orElseThrow(() -> new IllegalStateException("Cannot find " + rsid));

    String[] alleles = new String[2];
    for (int x = 0; x < calls.length; x += 1) {
      if (calls[x].equals("0")) {
        alleles[x] = vl.getRef();
      } else {
        alleles[x] = nonPgxAllele == null ? vl.getAlts().get(0) : nonPgxAllele;
      }
    }
    vcfBuilder.variation("DPYD", rsid, alleles[0], alleles[1]);
    if (nonPgxAllele != null) {
      vcfBuilder.allowUnknownAllele();
    }
  }


  /**
   * This test and {@link #hapB3AndIntronicB(TestInfo)} are to make sure we call diplotypes with HapB3 on one strand and
   * HapB3Intronic on the other correctly.
   * <p>
   * This test has GT calls on opposite strands from {@link #hapB3AndIntronicB(TestInfo)}.
   */
  @Test
  void hapB3AndIntronicA(TestInfo testInfo) throws Exception {

    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, true, false, false);
    testWrapper.getVcfBuilder()
        .phased()
        // hapB3 exon C>T
        .variation("DPYD", "rs56038477", "T", "C")
        // hapB3 intron G>C
        .variation("DPYD", "rs75017182", "C", "C")
        // *2A C>T
        .variation("DPYD", "rs3918290", "T", "C")
    ;

    Path vcfFile = testWrapper.execute();

    List<String> expectedCalls = List.of(
        "c.1129-5923C>G/[c.1129-5923C>G, c.1236G>A (HapB3) + c.1905+1G>A (*2A)]"
    );

    doStandardChecks(testWrapper, vcfFile, expectedCalls, null, List.of("c.1905+1G>A (*2A)", "c.1129-5923C>G"),
        false, RecPresence.YES);
  }

  /**
   * This test and {@link #hapB3AndIntronicA(TestInfo)} are to make sure we call diplotypes with HapB3 on one strand and
   * HapB3Intronic on the other correctly.
   * <p>
   * This test has GT calls on opposite strands from {@link #hapB3AndIntronicA(TestInfo)}.
   */
  @Test
  void hapB3AndIntronicB(TestInfo testInfo) throws Exception {

    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, true, false, false);
    testWrapper.getVcfBuilder()
        .phased()
        // hapB3 exon C>T
        .variation("DPYD", "rs56038477", "C", "T")
        // hapB3 intron G>C
        .variation("DPYD", "rs75017182", "C", "C")
        // *2A C>T
        .variation("DPYD", "rs3918290", "C", "T")
    ;

    Path vcfFile = testWrapper.execute();

    List<String> expectedCalls = List.of(
        "c.1129-5923C>G/[c.1129-5923C>G, c.1236G>A (HapB3) + c.1905+1G>A (*2A)]"
    );
    List<String> recommendedDiplotypes = List.of("c.1905+1G>A (*2A)", "c.1129-5923C>G");

    doStandardChecks(testWrapper, vcfFile, expectedCalls, null, recommendedDiplotypes, false, RecPresence.YES);
  }
}

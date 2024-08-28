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

  public static void dpydHasReports(PipelineWrapper testWrapper, RecPresence hasCpicReport, RecPresence hasDpwgReport) {
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



  @Test
  void testDpydUndocumentedVariation(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .allowUnknownAllele()
        .variation("DPYD", "rs3918290", "C", "G")
    ;
    Path vcfFile = testWrapper.execute();

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


  private static void dpydHtmlChecks(Document document, @Nullable List<String> expectedCalls,
      boolean hasMissingPositions, RecPresence hasDpwgAnnotation) {
    dpydHtmlChecks(document, expectedCalls, null, hasMissingPositions, RecPresence.YES, hasDpwgAnnotation);
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
        assertEquals(1, missingVariantsWarning.size());
      } else {
        assertEquals(0, missingVariantsWarning.size());
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
  void testDpydC2846het(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .phased()
        .variation("DPYD", "rs67376798", "T", "A");
    Path vcfFile = testWrapper.execute();

    List<String> expectedCalls = List.of("Reference/c.2846A>T");
    List<String> cpicStyleCalls = List.of("c.2846A>T (heterozygous)");
    RecPresence hasDpwgAnnotations = RecPresence.YES;

    testWrapper.testCalledByMatcher("DPYD");
    testWrapper.testSourceDiplotypes(DataSource.CPIC, "DPYD", expectedCalls, cpicStyleCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "DPYD", expectedCallsToRecommendedDiplotypes(expectedCalls));
    testWrapper.testPrintCalls(DataSource.CPIC, "DPYD", cpicStyleCalls);

    dpydHasReports(testWrapper, hasDpwgAnnotations);

    Document document = readHtmlReport(vcfFile);
    dpydHtmlChecks(document, expectedCalls, cpicStyleCalls, false, RecPresence.YES, hasDpwgAnnotations);
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


  @Test
  void testDpydS12het(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .phased()
        .variation("DPYD", "rs78060119", "C", "A");
    Path vcfFile = testWrapper.execute();

    List<String> expectedCalls = List.of("Reference/c.1156G>T (*12)");
    List<String> cpicStyleCalls = List.of("c.1156G>T (*12) (heterozygous)");
    RecPresence hasDpwgAnnotations = RecPresence.NO;

    testWrapper.testCalledByMatcher("DPYD");
    testWrapper.testSourceDiplotypes(DataSource.CPIC, "DPYD", cpicStyleCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "DPYD", expectedCallsToRecommendedDiplotypes(expectedCalls));
    testWrapper.testPrintCalls(DataSource.CPIC, "DPYD", cpicStyleCalls);

    dpydHasReports(testWrapper, hasDpwgAnnotations);

    Document document = readHtmlReport(vcfFile);
    dpydHtmlChecks(document, expectedCalls, cpicStyleCalls, false, RecPresence.YES, hasDpwgAnnotations);
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
    RecPresence hasDpwgAnnotations = RecPresence.YES;

    testWrapper.testCalledByMatcher("DPYD");
    testWrapper.testSourceDiplotypes(DataSource.CPIC, "DPYD", expectedCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "DPYD", List.of("c.1129-5923C>G, c.1236G>A (HapB3)", "c.1129-5923C>G"));
    testWrapper.testPrintCalls(DataSource.CPIC, "DPYD", expectedCalls);

    dpydHasReports(testWrapper, hasDpwgAnnotations);

    Document document = readHtmlReport(vcfFile);
    dpydHtmlChecks(document, expectedCalls, false, hasDpwgAnnotations);
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
    RecPresence hasDpwgAnnotations = RecPresence.YES;

    testWrapper.testCalledByMatcher("DPYD");
    testWrapper.testSourceDiplotypes(DataSource.CPIC, "DPYD", expectedCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "DPYD", List.of("c.1129-5923C>G, c.1236G>A (HapB3)", "c.1129-5923C>G"));
    testWrapper.testPrintCalls(DataSource.CPIC, "DPYD", expectedCalls);

    dpydHasReports(testWrapper, hasDpwgAnnotations);

    Document document = readHtmlReport(vcfFile);
    dpydHtmlChecks(document, expectedCalls, false, hasDpwgAnnotations);
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
      assertEquals("rs1801266", data[4]);
      assertEquals("rs186169810", data[5]);
      assertEquals("rs115232898", data[6]);
      int row = 3;
      while ((line = reader.readLine()) != null) {
        row += 1;
        if (row < 107) {
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
    setVariation(testWrapper, vcfBuilder, "rs1801266", data[4]);
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
    RecPresence hasDpwgAnnotations = RecPresence.YES;

    testWrapper.testCalledByMatcher("DPYD");
    testWrapper.testSourceDiplotypes(DataSource.CPIC, "DPYD", expectedCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "DPYD", List.of("c.1905+1G>A (*2A)", "c.1129-5923C>G"));
    testWrapper.testPrintCalls(DataSource.CPIC, "DPYD", expectedCalls);

    dpydHasReports(testWrapper, hasDpwgAnnotations);

    Document document = readHtmlReport(vcfFile);
    dpydHtmlChecks(document, expectedCalls, false, hasDpwgAnnotations);
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
    RecPresence hasDpwgAnnotations = RecPresence.YES;

    testWrapper.testCalledByMatcher("DPYD");
    testWrapper.testSourceDiplotypes(DataSource.CPIC, "DPYD", expectedCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "DPYD", List.of("c.1905+1G>A (*2A)", "c.1129-5923C>G"));
    testWrapper.testPrintCalls(DataSource.CPIC, "DPYD", expectedCalls);

    dpydHasReports(testWrapper, hasDpwgAnnotations);

    Document document = readHtmlReport(vcfFile);
    dpydHtmlChecks(document, expectedCalls, false, hasDpwgAnnotations);
  }
}

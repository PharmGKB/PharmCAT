package org.pharmgkb.pharmcat;

import java.io.BufferedWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.stream.Collectors;
import com.google.common.base.Preconditions;
import com.google.common.collect.ImmutableSortedMap;
import com.google.common.collect.Ordering;
import org.checkerframework.checker.nullness.qual.Nullable;
import org.jsoup.Jsoup;
import org.jsoup.nodes.Document;
import org.jsoup.nodes.Element;
import org.jsoup.select.Elements;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.TestInfo;
import org.pharmgkb.common.util.PathUtils;
import org.pharmgkb.pharmcat.reporter.MessageHelper;
import org.pharmgkb.pharmcat.reporter.TextConstants;
import org.pharmgkb.pharmcat.reporter.handlebars.ReportHelpers;
import org.pharmgkb.pharmcat.reporter.model.DataSource;
import org.pharmgkb.pharmcat.reporter.model.MessageAnnotation;
import org.pharmgkb.pharmcat.reporter.model.PrescribingGuidanceSource;
import org.pharmgkb.pharmcat.reporter.model.VariantReport;
import org.pharmgkb.pharmcat.reporter.model.result.AnnotationReport;
import org.pharmgkb.pharmcat.reporter.model.result.Diplotype;
import org.pharmgkb.pharmcat.reporter.model.result.DrugReport;
import org.pharmgkb.pharmcat.reporter.model.result.GeneReport;
import org.pharmgkb.pharmcat.reporter.model.result.GuidelineReport;
import org.pharmgkb.pharmcat.reporter.model.result.Haplotype;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.contains;
import static org.junit.jupiter.api.Assertions.*;
import static org.pharmgkb.pharmcat.Constants.isLowestFunctionGene;


/**
 * This is a JUnit test for {@link Pipeline}.
 * This should test the data generated from a full run of the PharmCAT matcher and reporter.
 *
 * @author Mark Woon
 */
class PipelineTest {
  private static final String sf_unknownDiplotype = Haplotype.UNKNOWN + TextConstants.GENOTYPE_DELIMITER + Haplotype.UNKNOWN;
  public static final List<String> UNKNOWN_CALL = List.of(sf_unknownDiplotype);
  public static final List<String> NO_DATA = List.of(sf_unknownDiplotype);
  public static final List<String> NO_OUTSIDE_DIPLOTYPE = List.of(sf_unknownDiplotype);
  private static Path s_outsideCallFilePath;
  private static Path s_otherOutsideCallFilePath;


  @BeforeAll
  static void prepare() throws IOException {

    ReportHelpers.setDebugMode(true);
    s_outsideCallFilePath = TestUtils.createTestFile(PipelineTest.class, "outsideCall.tsv");
    try (BufferedWriter writer = Files.newBufferedWriter(s_outsideCallFilePath)) {
      writer.write("""
          ##Test Outside Call Data
          CYP2D6\t*1/*4\t\t\t0.6\t0.75\tp: 0.0\t\t\tv1.9-2017_02_09
          """);
    }

    s_otherOutsideCallFilePath = TestUtils.createTestFile(PipelineTest.class, "otherOutsideCall.tsv");
    try (BufferedWriter writer = Files.newBufferedWriter(s_otherOutsideCallFilePath)) {
      writer.write("""
          CYP2D6\t*3/*4
          G6PD\tB (wildtype)/B (wildtype)
          """);
    }
    TestUtils.setSaveTestOutput(true);
  }

  @AfterEach
  void deleteDirectory(TestInfo testInfo) {
    TestUtils.deleteTestOutputDirectory(testInfo);
  }


  public static List<String> expectedCallsToRecommendedDiplotypes(List<String> expectedCalls) {
    Preconditions.checkArgument(expectedCalls.size() == 1);
    return expectedCalls.stream()
        .flatMap(s -> Arrays.stream(s.split(TextConstants.GENOTYPE_DELIMITER)))
        .toList();
  }


  /**
   * Checks for expected HTML output.
   *
   * @param expectedCalls - use {@link #UNKNOWN_CALL} or {@link #NO_DATA} where necessary
   */
  public static void htmlChecks(Document document, String gene, List<String> expectedCalls,
      String drug, RecPresence cpicAnnPresence, RecPresence dpwgAnnPresence) {
    Preconditions.checkNotNull(expectedCalls);
    htmlChecks(document,
        new ImmutableSortedMap.Builder<String, List<String>>(Ordering.natural())
            .put(gene, expectedCalls)
            .build(),
        null, drug, cpicAnnPresence, dpwgAnnPresence);
  }

  /**
   * Checks for expected HTML output.
   *
   * @param expectedCalls - use {@link #UNKNOWN_CALL} or {@link #NO_DATA} where necessary
   */
  public static void htmlChecks(Document document, SortedMap<String, List<String>> expectedCalls,
      String drug, RecPresence cpicAnnPresence, RecPresence dpwgAnnPresence) {
    htmlChecks(document, expectedCalls, null, drug, cpicAnnPresence, dpwgAnnPresence);
  }

  /**
   * Checks for expected HTML output.
   *
   * @param expectedCalls - use {@link #UNKNOWN_CALL} or {@link #NO_DATA} where necessary
   */
  static void htmlChecks(Document document, SortedMap<String, List<String>> expectedCalls,
      @Nullable Map<String, List<String>> cpicStyleCalls, String drug, RecPresence cpicAnnPresence,
      RecPresence dpwgAnnPresence) {
    htmlChecks(document, expectedCalls, cpicStyleCalls, drug, cpicAnnPresence, null, dpwgAnnPresence, null);
  }

  /**
   * Checks for expected HTML output.
   *
   * @param expectedCalls - use {@link #UNKNOWN_CALL} or {@link #NO_DATA} where necessary
   */
  static void htmlChecks(Document document, SortedMap<String, List<String>> expectedCalls,
      @Nullable Map<String, List<String>> cpicStyleCalls, String drug, RecPresence cpicAnnPresence,
      @Nullable SortedMap<String, String> cpicPhenotypes, RecPresence dpwgAnnPresence,
      @Nullable SortedMap<String, String> dpwgPhenotypes) {
    htmlChecks(document, expectedCalls, cpicStyleCalls, drug, cpicAnnPresence, cpicPhenotypes, null, dpwgAnnPresence,
        dpwgPhenotypes, null);
  }

  /**
   * Checks for expected HTML output.
   *
   * @param expectedCalls - use {@link #UNKNOWN_CALL} or {@link #NO_DATA} where necessary
   */
  static void htmlChecks(Document document, SortedMap<String, List<String>> expectedCalls,
      @Nullable Map<String, List<String>> cpicStyleCalls, String drug,
      RecPresence cpicAnnPresence,
      @Nullable SortedMap<String, String> cpicPhenotypes,
      @Nullable SortedMap<String, SortedSet<String>> cpicActivityScores,
      RecPresence dpwgAnnPresence,
      @Nullable SortedMap<String, String> dpwgPhenotypes,
      @Nullable SortedMap<String, SortedSet<String>> dpwgActivityScores) {
    if (drug == null) {
      throw new IllegalArgumentException("If you don't want to provide drug, use htmlCheckGenes()");
    }

    htmlCheckGenes(document, expectedCalls, cpicStyleCalls);
    htmlCheckDrug(document, expectedCalls, cpicStyleCalls, drug, cpicAnnPresence, cpicPhenotypes, cpicActivityScores,
        dpwgAnnPresence, dpwgPhenotypes, dpwgActivityScores);
  }

  /**
   * Runs all gene-related HTML checks.
   *
   * @param expectedCalls - use {@link #UNKNOWN_CALL} or {@link #NO_DATA} where necessary
   */
  static void htmlCheckGenes(Document document, Map<String, List<String>> expectedCalls,
      @Nullable Map<String, List<String>> cpicStyleCalls) {
    for (String gene : expectedCalls.keySet()) {
      htmlCheckGene(document, gene, expectedCalls.get(gene), cpicStyleCalls == null ? null : cpicStyleCalls.get(gene));
    }
  }

  static void htmlCheckGene(Document document, String gene, List<String> expectedCalls,
      @Nullable List<String> cpicStyleCalls) {
    Preconditions.checkNotNull(expectedCalls);
    if (cpicStyleCalls != null && cpicStyleCalls.isEmpty()) {
      cpicStyleCalls = null;
    }

    if (expectedCalls == NO_DATA || expectedCalls == UNKNOWN_CALL) {
      // check section i
      assertEquals(0, document.select(".gs-" + gene + " .gs-dip").size());
      if (expectedCalls == UNKNOWN_CALL) {
        assertNotNull(document.getElementById("gs-uncallable-" + gene), gene + " should be uncallable");
      }

      // check section iii
      Elements geneSection = document.select(".gene." + gene);
      assertEquals(1, geneSection.size());
      if (expectedCalls == NO_DATA) {
        assertEquals(1, geneSection.get(0).getElementsByClass("no-data").size());
      } else {
        assertEquals("Not called", geneSection.select(".genotype-result").text());
      }

    } else {
      // check section i
      boolean didLowestFunctionCheck = false;
      if (isLowestFunctionGene(gene)) {
        List<String> expectedComponents = DpydTest.callsToComponents(expectedCalls);

        if (expectedComponents != null) {
          Elements gsLowestFunction = document.select(".gs-" + gene + " .gs-dip_lowestFunction");
          assertEquals(cpicStyleCalls == null ? expectedCalls : cpicStyleCalls,
              gsLowestFunction.stream()
                  .map(e -> e.child(0).text())
                  .toList());

          Elements gsComponents = document.select(".gs-" + gene + " .gs-dip_component");
          List<String> components = gsComponents.stream()
              .map(e -> e.child(0).text())
              .toList();
          assertEquals(expectedComponents, components);
          didLowestFunctionCheck = true;
        }
      }
      if (!didLowestFunctionCheck) {
        List<String> gsDips = document.select(".gs-" + gene + " .gs-dip").stream()
            .map(e -> e.child(0).text())
            .toList();
        if (expectedCalls == NO_OUTSIDE_DIPLOTYPE) {
          assertEquals(1, gsDips.size());
          assertEquals(TextConstants.OUTSIDE_DATA_NO_GENOTYPE, gsDips.get(0));
        } else {
          assertEquals(cpicStyleCalls == null ? expectedCalls : cpicStyleCalls, gsDips);
        }
      }
      // check section iii
      Elements geneSection = document.select(".gene." + gene);
      assertEquals(1, geneSection.size());
      assertEquals(0, geneSection.get(0).getElementsByClass("no-data").size());
    }
  }

  private static void htmlCheckDrug(Document document, SortedMap<String, List<String>> expectedCalls,
      @Nullable Map<String, List<String>> cpicStyleCalls, String drug, RecPresence cpicAnnPresence,
      @Nullable SortedMap<String, String> cpicPhenotypes, @Nullable SortedMap<String, SortedSet<String>> cpicActivityScores,
      RecPresence dpwgAnnPresence, @Nullable SortedMap<String, String> dpwgPhenotypes,
      @Nullable SortedMap<String, SortedSet<String>> dpwgActivityScores) {

    String sanitizedDrug = ReportHelpers.sanitizeCssSelector(drug);
    Elements drugSections = document.getElementsByClass(sanitizedDrug);

    if (cpicAnnPresence == RecPresence.NO && dpwgAnnPresence == RecPresence.NO) {
      assertEquals(0, drugSections.size());

    } else {
      assertEquals(1, drugSections.size());

      SortedMap<String, List<String>> expectedRxCalls = new TreeMap<>();
      for (String gene : expectedCalls.keySet()) {
        List<String> calls = expectedCalls.get(gene);
        if (cpicStyleCalls != null && cpicStyleCalls.containsKey(gene)) {
          calls = cpicStyleCalls.get(gene);
        }
        if (calls == null || calls.isEmpty()) {
          Elements cpicDrugDips = drugSections.select(".cpic-guideline-" + sanitizedDrug + " .rx-dip");
          assertEquals(0, cpicDrugDips.size());

          Elements dpwgDrugDips = drugSections.select(".dpwg-guideline-" + sanitizedDrug + " .rx-dip");
          assertEquals(0, dpwgDrugDips.size());

          continue;
        }

        if (calls == NO_DATA) {
          expectedRxCalls.put(gene, List.of(gene + ":" + TextConstants.NO_DATA));
        } else if (calls == NO_OUTSIDE_DIPLOTYPE) {
          expectedRxCalls.put(gene, List.of(gene + ":" + TextConstants.OUTSIDE_DATA_NO_GENOTYPE));
        } else {
          expectedRxCalls.put(gene, calls.stream()
              .map(c -> gene + ":" + c)
              .collect(Collectors.toList()));
        }
      }

      htmlCheckDrugAnnotation(drugSections, "cpic-guideline", sanitizedDrug, cpicAnnPresence,
          filterRxCalls(expectedRxCalls, cpicPhenotypes), cpicPhenotypes, cpicActivityScores);

      htmlCheckDrugAnnotation(drugSections, "dpwg-guideline", sanitizedDrug, dpwgAnnPresence,
          filterRxCalls(expectedRxCalls, dpwgPhenotypes), dpwgPhenotypes, dpwgActivityScores);
    }
  }

  private static SortedMap<String, List<String>> filterRxCalls(SortedMap<String, List<String>> expectedRxCalls,
      Map<String, String> phenotypes) {
    if (phenotypes == null || phenotypes.size() == expectedRxCalls.size()) {
      return expectedRxCalls;
    }
    SortedMap<String, List<String>> filteredRxCalls = new TreeMap<>();
    for (String gene : phenotypes.keySet()) {
      filteredRxCalls.put(gene, expectedRxCalls.get(gene));
    }
    return filteredRxCalls;
  }


  private static void htmlCheckDrugAnnotation(Elements drugSections, String src, String drug, RecPresence annPresence,
      SortedMap<String, List<String>> expectedRxCalls, @Nullable SortedMap<String, String> expectedPhenotypes,
      @Nullable SortedMap<String, SortedSet<String>> expectedActivityScores) {

    String srcSelector = "." + src + "-" + ReportHelpers.sanitizeCssSelector(drug);
    Elements srcSections = drugSections.select(srcSelector);

    if (annPresence == RecPresence.NO) {
      Elements rows = srcSections.select(srcSelector + " .rx-dip");
      assertEquals(0, rows.size());
      rows = srcSections.select(srcSelector + " .rx-unmatched-dip");
      assertEquals(0, rows.size());

    } else if (annPresence == RecPresence.YES_NO_MATCH) {
      Elements unmatchedDips = drugSections.select(srcSelector + " .rx-unmatched-dip");
      assertEquals(expectedRxCalls.values().stream()
              .flatMap(Collection::stream)
              .toList(),
          unmatchedDips.stream()
              .map(e -> cleanupRxDip(e, expectedRxCalls.keySet()))
              .toList());

    } else {
      SortedMap<String, SortedSet<String>> actualActivityScores = new TreeMap<>();
      boolean hasActivityScore = false;
      for (Element row : srcSections) {
        Elements rxDips = row.select(".rx-dip");
        assertEquals(expectedRxCalls.values().stream()
                .flatMap(Collection::stream)
                .toList(),
            rxDips.stream()
                .map(e -> cleanupRxDip(e, expectedRxCalls.keySet()))
                .toList());

        Elements rxPhenotypes = row.select(".rx-phenotype");
        if (expectedPhenotypes != null) {
          assertEquals(expectedPhenotypes.size(), rxPhenotypes.size());
          SortedMap<String, String> actualPhenotypes = new TreeMap<>();
          if (rxPhenotypes.size() == 1) {
            actualPhenotypes.put(expectedPhenotypes.firstKey(), rxPhenotypes.get(0).text());
          } else {
            for (Element e : rxPhenotypes) {
              Elements dts = e.select("dt");
              String gene = dts.get(0).text();
              if (gene.endsWith(":")) {
                gene = gene.substring(0, gene.length() - 1);
              }
              actualPhenotypes.put(gene, e.select("dd").get(0).text());
            }
          }
          assertEquals(expectedPhenotypes, actualPhenotypes);
        }

        Elements rxActivityScores = row.select(".rx-activity");
        if (!rxActivityScores.isEmpty()) {
          hasActivityScore = true;
        }
        if (expectedActivityScores != null) {
          assertEquals(expectedActivityScores.size(), rxActivityScores.size());
          if (rxActivityScores.size() == 1) {
            String gene = expectedActivityScores.firstKey();
            actualActivityScores.putIfAbsent(gene, new TreeSet<>());
            actualActivityScores.get(gene).add(rxActivityScores.get(0).text());
          } else {
            for (Element e : rxActivityScores) {
              Elements dts = e.select("dt");
              String gene = dts.get(0).text();
              if (gene.endsWith(":")) {
                gene = gene.substring(0, gene.length() - 1);
              }
              actualActivityScores.putIfAbsent(gene, new TreeSet<>());
              actualActivityScores.get(gene).add(e.select("dd").get(0).text());
            }
          }
        }
      }

      // we check activity scores across all rows
      if (expectedActivityScores == null) {
        if (hasActivityScore) {
          fail("Not checking for activity score!");
        }
      } else {
        assertEquals(expectedActivityScores, actualActivityScores);
      }
    }
  }

  public static String cleanupRxDip(Element rxDip, Collection<String> genes) {
    String dip = rxDip.text().replace("/ ", "/");
    for (String gene : genes) {
      dip = dip.replace(gene + ": ", gene + ":");
    }
    return dip;
  }


  /**
   * NOTE: if these assertions fail, then new data may have been added from the DataManager because of an update to the
   * CPIC database. If that's true, then update these numbers to the current count. If the count changes with no known
   * change to the CPIC database, then something may be wrong in the code.
   *
   * <p>NOTE: you may also want to check the git log of Genes-Drugs.md to see specifically what genes/drugs have
   * changed.</p>
   */
  @Test
  void testCounts(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("CYP2C9");
    testWrapper.execute();
    SortedSet<String> genes = testWrapper.getContext().getGeneReports().keySet().stream()
        .flatMap((k) -> testWrapper.getContext().getGeneReports().get(k).values().stream()
            .map(GeneReport::getGeneDisplay))
        .collect(Collectors.toCollection(TreeSet::new));
    assertEquals(22, genes.size());

    SortedSet<String> drugs = testWrapper.getContext().getDrugReports().keySet().stream()
        .flatMap((k) -> testWrapper.getContext().getDrugReports().get(k).values().stream()
            .map(DrugReport::getName))
        .collect(Collectors.toCollection(TreeSet::new));
    assertEquals(183, drugs.size());
  }

  @Test
  void testAll(TestInfo testInfo) throws Exception {
    Path outsideCallPath = TestUtils.createTestFile(testInfo, ".tsv");
    try (PrintWriter writer = new PrintWriter(Files.newBufferedWriter(outsideCallPath))) {
      writer.println(
          """
              CYP2D6\t*3/*4
              HLA-A\t\t*31:01 positive
              HLA-B\t*15:02/*57:01"""
      );
    }
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("ABCG2")
        .reference("CACNA1S")
        .reference("CFTR")
        .reference("CYP2B6")
        .reference("CYP2C19")
        .variation("CYP2C19", "rs3758581", "G", "G") // to make it *1/*1
        .reference("CYP2C9")
        .reference("CYP3A4")
        .reference("CYP3A5")
        .reference("CYP4F2")
        .reference("DPYD")
        .reference("G6PD")
        .reference("IFNL3")
        .reference("NUDT15")
        .reference("RYR1")
        .reference("SLCO1B1")
        .reference("TPMT")
        .reference("UGT1A1")
        .reference("VKORC1");
    testWrapper.executeWithOutsideCalls(outsideCallPath);

    testWrapper.testCalledByMatcher(
        "ABCG2",
        "CACNA1S",
        "CFTR",
        "CYP2B6",
        "CYP2C19",
        "CYP2C9",
        "CYP3A4",
        "CYP3A5",
        "CYP4F2",
        "DPYD",
        "G6PD",
        "IFNL3",
        "NUDT15",
        "RYR1",
        "SLCO1B1",
        "TPMT",
        "UGT1A1",
        "VKORC1"
    );
    testWrapper.testNotCalledByMatcher("CYP2D6", "HLA-A", "HLA-B");
  }


  @Test
  void testNoData(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.execute(null,  null, null, true);
  }


  @Test
  void testUndocumentedVariation(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .allowUnknownAllele()
        .variation("CYP2C19", "rs3758581", "G", "T");
    Path vcfFile = testWrapper.execute();

    testWrapper.testNotCalledByMatcher("CYP2C19");

    Document document = readHtmlReport(vcfFile);
    htmlCheckGenes(document,
        new ImmutableSortedMap.Builder<String, List<String>>(Ordering.natural())
            .put("CYP2C19", UNKNOWN_CALL)
            .build(),
        null);
  }

  @Test
  void testUndocumentedVariationExtendedReport(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false)
        .extendedReport();
    testWrapper.getVcfBuilder()
        .allowUnknownAllele()
        .variation("CYP2C19", "rs3758581", "G", "T");
    Path vcfFile = testWrapper.execute();

    testWrapper.testNotCalledByMatcher("CYP2C19");

    Document document = readHtmlReport(vcfFile);
    htmlCheckGenes(document,
        new ImmutableSortedMap.Builder<String, List<String>>(Ordering.natural())
            .put("CYP2C19", UNKNOWN_CALL)
            .build(),
        null);
  }

  @Test
  void testUndocumentedVariationsWithTreatAsReference(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .allowUnknownAllele()
        // undocumented as reference
        .variation("TPMT", "rs1800462", "C", "T")
        // undocumented as reference + lowest-function
        .variation("RYR1", "rs193922753", "G", "C")
        // not undocumented-as-reference
        .variation("CYP2C19", "rs3758581", "G", "T");
    Path vcfFile = testWrapper.execute();

    List<String> tpmtExpectedCalls = List.of("*1/*1");
    List<String> ryr1ExpectedCalls = List.of(TextConstants.HOMOZYGOUS_REFERENCE);

    testWrapper.testNotCalledByMatcher("CYP2C19");
    testWrapper.testCalledByMatcher("TPMT");
    testWrapper.testCalledByMatcher("RYR1");
    testWrapper.testSourceDiplotypes(DataSource.CPIC, "TPMT", tpmtExpectedCalls);
    testWrapper.testSourceDiplotypes(DataSource.CPIC, "RYR1", ryr1ExpectedCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "TPMT", expectedCallsToRecommendedDiplotypes(tpmtExpectedCalls));
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "RYR1", expectedCallsToRecommendedDiplotypes(ryr1ExpectedCalls));

    Document document = readHtmlReport(vcfFile);
    assertNotNull(document.getElementById("gs-undocVarAsRef-TPMT"));
    assertNotNull(document.getElementById("gs-undocVarAsRef-RYR1"));

    htmlCheckGenes(document,
        new ImmutableSortedMap.Builder<String, List<String>>(Ordering.natural())
            .put("CYP2C19", UNKNOWN_CALL)
            .put("TPMT", tpmtExpectedCalls)
            .put("RYR1", ryr1ExpectedCalls)
            .build(),
        null);
  }

  @Test
  void testUndocumentedVariationsWithTreatAsReferenceAndCombo(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, true, false, false);
    testWrapper.getVcfBuilder()
        .allowUnknownAllele()
        // undocumented as reference
        .variation("TPMT", "rs1800462", "C", "T")
        // undocumented as reference + lowest-function
        .variation("RYR1", "rs193922753", "G", "C");
    Path vcfFile = testWrapper.execute();

    testWrapper.testCalledByMatcher("TPMT");
    testWrapper.testCalledByMatcher("RYR1");

    // becomes Reference and custom snp because combo is enabled
    List<String> tpmtExpectedCalls = List.of("*1/g.18143724C>T");
    // but lowest-function ignores combo
    List<String> ryr1ExpectedCalls = List.of(TextConstants.HOMOZYGOUS_REFERENCE);

    testWrapper.testSourceDiplotypes(DataSource.CPIC, "TPMT", tpmtExpectedCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "TPMT", expectedCallsToRecommendedDiplotypes(tpmtExpectedCalls));
    testWrapper.testSourceDiplotypes(DataSource.CPIC, "RYR1", ryr1ExpectedCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "RYR1", expectedCallsToRecommendedDiplotypes(ryr1ExpectedCalls));
    testWrapper.testRecommendedPhenotype(DataSource.CPIC, "TPMT", TextConstants.NA);

    Document document = readHtmlReport(vcfFile);
    assertNull(document.getElementById("gs-undocVarAsRef-TPMT"));
    assertNotNull(document.getElementById("gs-undocVarAsRef-RYR1"));

    htmlCheckGenes(document,
        new ImmutableSortedMap.Builder<String, List<String>>(Ordering.natural())
            .put("RYR1", ryr1ExpectedCalls)
            .build(),
        null);
  }


  @Test
  void testUncallable(TestInfo testInfo) throws Exception {

    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("ABCG2")
        .allowUnknownAllele()
        .variation("CYP2C19", "rs3758581", "G", "T")
        .variation("TPMT", "rs1256618794", "A", "A") // C -> A
        .variation("TPMT", "rs753545734", "C", "C") // C -> T
    ;
    Path vcfFile = testWrapper.execute();

    List<String> expectedCalls = UNKNOWN_CALL;

    testWrapper.testNotCalledByMatcher("CYP2C19", "TPMT");

    testWrapper.testSourceDiplotypes(DataSource.CPIC, "CYP2C19", expectedCalls);
    testWrapper.testSourceDiplotypes(DataSource.CPIC, "TPMT", expectedCalls);

    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "CYP2C19", expectedCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "TPMT", expectedCalls);

    testWrapper.testPrintCalls(DataSource.CPIC, "CYP2C19", List.of(TextConstants.UNCALLED));
    testWrapper.testPrintCalls(DataSource.CPIC, "TPMT", List.of(TextConstants.UNCALLED));

    Document document = readHtmlReport(vcfFile);

    SortedMap<String, List<String>> expectedCallsMap = new TreeMap<>();
    expectedCallsMap.put("CYP2C19", UNKNOWN_CALL);
    expectedCallsMap.put("TPMT", UNKNOWN_CALL);
    htmlCheckGenes(document, expectedCallsMap, null);
  }


  /**
   * This test illustrates when one gene in a two-gene guideline (amitriptyline) is not called that it should still be
   * able to come up with a matched annotation.
   */
  @Test
  void testCyp2c19(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("CYP2C19", "rs3758581", "G", "G");
    testWrapper.executeWithOutsideCalls(s_otherOutsideCallFilePath);

    testWrapper.testCalledByMatcher("CYP2C19");
    testWrapper.testPrintCpicCalls( "CYP2C19", "*1/*1");

    testWrapper.testMatchedAnnotations("amitriptyline", PrescribingGuidanceSource.CPIC_GUIDELINE, 1);
    testWrapper.testMatchedAnnotations("amitriptyline", PrescribingGuidanceSource.DPWG_GUIDELINE, 1);
    testWrapper.testMatchedAnnotations("citalopram", PrescribingGuidanceSource.CPIC_GUIDELINE, 1);
    testWrapper.testMatchedAnnotations("citalopram", PrescribingGuidanceSource.DPWG_GUIDELINE, 1);
    testWrapper.testMatchedAnnotations("ivacaftor", 0);
  }

  /**
   * This test case demos that an "ambiguity" {@link MessageAnnotation} which specifies a variant and a diplotype call
   * for a given drug report will be matched and added to the {@link DrugReport}
   */
  @Test
  void testCyp2c19_s1s2rs58973490het(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("CYP2C19", "rs12769205", "A", "G")
        .variation("CYP2C19", "rs58973490", "G", "A")
        .variation("CYP2C19", "rs4244285", "G", "A")
        .variation("CYP2C19", "rs3758581", "G", "G");
    testWrapper.executeWithOutsideCalls(s_otherOutsideCallFilePath);

    testWrapper.testCalledByMatcher("CYP2C19");
    testWrapper.testPrintCpicCalls( "CYP2C19", "*1/*2");
    testWrapper.testNotCalledByMatcher("CYP2D6");
    testWrapper.testPrintCpicCalls( "CYP2D6", "*3/*4");

    testWrapper.testMatchedAnnotations("amitriptyline", 3);
    testWrapper.testAnyMatchFromSource("amitriptyline", PrescribingGuidanceSource.CPIC_GUIDELINE);
    testWrapper.testAnyMatchFromSource("amitriptyline", PrescribingGuidanceSource.DPWG_GUIDELINE);
    testWrapper.testAnyMatchFromSource("amitriptyline", PrescribingGuidanceSource.FDA_ASSOC);
    testWrapper.testMatchedAnnotations("citalopram", 2);
    testWrapper.testAnyMatchFromSource("citalopram", PrescribingGuidanceSource.CPIC_GUIDELINE);
    testWrapper.testAnyMatchFromSource("citalopram", PrescribingGuidanceSource.DPWG_GUIDELINE);
    testWrapper.testMatchedAnnotations("clomipramine", 4);
    testWrapper.testMatchedAnnotations("clomipramine", PrescribingGuidanceSource.CPIC_GUIDELINE, 1);
    testWrapper.testMatchedAnnotations("clomipramine", PrescribingGuidanceSource.DPWG_GUIDELINE, 2);
    testWrapper.testMatchedAnnotations("clomipramine", PrescribingGuidanceSource.FDA_ASSOC, 1);
    testWrapper.testMatchedAnnotations("ivacaftor", 0);

    GeneReport cyp2c19report = testWrapper.getContext().getGeneReport(DataSource.CPIC, "CYP2C19");
    assertNotNull(cyp2c19report);
    VariantReport vr = cyp2c19report.findVariantReport("rs58973490")
        .orElseThrow(() -> new RuntimeException("Variant missing from test data"));
    assertTrue(vr.isHetCall());

    // ambiguity message will not apply in this case because all variants are available for CYP2C19, but one message
    // should appear for the *1 call
    assertEquals(1, cyp2c19report.getMessages().stream()
        .filter(m -> m.getExceptionType().equals(MessageAnnotation.TYPE_AMBIGUITY) &&
            Objects.requireNonNull(m.getMatches().getVariant()).equals("rs58973490"))
        .count());

    testWrapper.testMessageCountForDrug(PrescribingGuidanceSource.CPIC_GUIDELINE, "amitriptyline", 1);
  }

  /**
   * This test case demos that an "ambiguity" {@link MessageAnnotation} which specifies a variant and a diplotype call
   * for a given drug report will not be matched when the variant in the message is homozygous
   */
  @Test
  void testCyp2c19_s1s2(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("CYP2C19", "rs12769205", "A", "G")
        .variation("CYP2C19", "rs4244285", "G", "A")
        .variation("CYP2C19", "rs3758581", "G", "G");
    testWrapper.executeWithOutsideCalls(s_otherOutsideCallFilePath);

    testWrapper.testCalledByMatcher("CYP2C19");
    testWrapper.testPrintCpicCalls( "CYP2C19", "*1/*2");

    testWrapper.testMatchedAnnotations("amitriptyline", 3);
    testWrapper.testAnyMatchFromSource("amitriptyline", PrescribingGuidanceSource.CPIC_GUIDELINE);
    testWrapper.testAnyMatchFromSource("amitriptyline", PrescribingGuidanceSource.DPWG_GUIDELINE);
    testWrapper.testAnyMatchFromSource("amitriptyline", PrescribingGuidanceSource.FDA_ASSOC);
    testWrapper.testMatchedAnnotations("citalopram", 2);
    testWrapper.testAnyMatchFromSource("citalopram", PrescribingGuidanceSource.CPIC_GUIDELINE);
    testWrapper.testAnyMatchFromSource("citalopram", PrescribingGuidanceSource.DPWG_GUIDELINE);
    testWrapper.testMatchedAnnotations("clomipramine", 4);
    testWrapper.testAnyMatchFromSource("clomipramine", PrescribingGuidanceSource.CPIC_GUIDELINE);
    testWrapper.testAnyMatchFromSource("clomipramine", PrescribingGuidanceSource.DPWG_GUIDELINE);
    testWrapper.testAnyMatchFromSource("clomipramine", PrescribingGuidanceSource.FDA_ASSOC);
    testWrapper.testMatchedAnnotations("ivacaftor", 0);

    GeneReport cyp2c19report = testWrapper.getContext().getGeneReport(DataSource.CPIC, "CYP2C19");
    assertNotNull(cyp2c19report);

    // make sure the variant in question is not a het call
    VariantReport vr = cyp2c19report.findVariantReport("rs58973490")
        .orElseThrow(() -> new RuntimeException("Variant missing from test data"));
    assertFalse(vr.isHetCall());

    // the variant is hom, so the ambiguity message should not apply and, thus, no matching messages
    assertEquals(0, cyp2c19report.getMessages().stream()
        .filter(m -> m.getExceptionType().equals(MessageAnnotation.TYPE_AMBIGUITY) &&
            Objects.requireNonNull(m.getMatches().getVariant()).equals("rs58973490"))
        .count());

    testWrapper.testAnyMatchFromSource("amitriptyline", PrescribingGuidanceSource.CPIC_GUIDELINE);
    // the variant is hom, so the ambiguity message should not apply and, thus, no matching messages
    testWrapper.testMessageCountForDrug(PrescribingGuidanceSource.CPIC_GUIDELINE, "amitriptyline", 0);

    // CYP2C19 reference is *38, not *1, so should not have a reference message
    testWrapper.testMessageCountForGene(DataSource.CPIC, "CYP2C19", 0);
  }

  @Test
  void testClomipramineCall(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("CYP2C19", "rs12769205", "G", "G")
        .variation("CYP2C19", "rs4244285", "A", "A")
        .variation("CYP2C19", "rs3758581", "G", "G");
    testWrapper.executeWithOutsideCalls(s_otherOutsideCallFilePath);

    testWrapper.testCalledByMatcher("CYP2C19");
    testWrapper.testPrintCpicCalls( "CYP2C19", "*2/*2");

    testWrapper.testMatchedAnnotations("amitriptyline", 3);
    testWrapper.testAnyMatchFromSource("amitriptyline", PrescribingGuidanceSource.CPIC_GUIDELINE);
    testWrapper.testAnyMatchFromSource("amitriptyline", PrescribingGuidanceSource.DPWG_GUIDELINE);
    testWrapper.testAnyMatchFromSource("amitriptyline", PrescribingGuidanceSource.FDA_ASSOC);
    testWrapper.testMatchedAnnotations("clomipramine", 4);
    testWrapper.testAnyMatchFromSource("clomipramine", PrescribingGuidanceSource.CPIC_GUIDELINE);
    testWrapper.testAnyMatchFromSource("clomipramine", PrescribingGuidanceSource.DPWG_GUIDELINE);
    testWrapper.testMatchedAnnotations("desipramine", 2);
    testWrapper.testAnyMatchFromSource("desipramine", PrescribingGuidanceSource.CPIC_GUIDELINE);
    testWrapper.testAnyMatchFromSource("desipramine", PrescribingGuidanceSource.FDA_ASSOC);
    testWrapper.testMatchedAnnotations("doxepin", 4);
    testWrapper.testAnyMatchFromSource("doxepin", PrescribingGuidanceSource.CPIC_GUIDELINE);
    testWrapper.testAnyMatchFromSource("doxepin", PrescribingGuidanceSource.DPWG_GUIDELINE);
    testWrapper.testAnyMatchFromSource("doxepin", PrescribingGuidanceSource.FDA_ASSOC);
    testWrapper.testMatchedAnnotations("imipramine", 4);
    testWrapper.testAnyMatchFromSource("imipramine", PrescribingGuidanceSource.CPIC_GUIDELINE);
    testWrapper.testAnyMatchFromSource("imipramine", PrescribingGuidanceSource.DPWG_GUIDELINE);
    testWrapper.testAnyMatchFromSource("imipramine", PrescribingGuidanceSource.FDA_ASSOC);
    testWrapper.testMatchedAnnotations("nortriptyline", 3);
    testWrapper.testAnyMatchFromSource("nortriptyline", PrescribingGuidanceSource.CPIC_GUIDELINE);
    testWrapper.testAnyMatchFromSource("nortriptyline", PrescribingGuidanceSource.DPWG_GUIDELINE);
    testWrapper.testAnyMatchFromSource("nortriptyline", PrescribingGuidanceSource.FDA_ASSOC);
    testWrapper.testMatchedAnnotations("trimipramine", 2);

    testWrapper.testMatchedAnnotations("clopidogrel", 6);
    testWrapper.testAnyMatchFromSource("clopidogrel", PrescribingGuidanceSource.CPIC_GUIDELINE);
    testWrapper.testAnyMatchFromSource("clopidogrel", PrescribingGuidanceSource.DPWG_GUIDELINE);
    testWrapper.testAnyMatchFromSource("clopidogrel", PrescribingGuidanceSource.FDA_LABEL);
    testWrapper.testAnyMatchFromSource("clopidogrel", PrescribingGuidanceSource.FDA_ASSOC);

    testWrapper.testMatchedAnnotations("lansoprazole", 3);
    testWrapper.testAnyMatchFromSource("lansoprazole", PrescribingGuidanceSource.CPIC_GUIDELINE);
    testWrapper.testAnyMatchFromSource("lansoprazole", PrescribingGuidanceSource.DPWG_GUIDELINE);
    testWrapper.testAnyMatchFromSource("lansoprazole", PrescribingGuidanceSource.FDA_ASSOC);

    // voriconazole has 2 populations with recommendations so should have 2 matching annotations from CPIC
    // and 1 from DPWG
    testWrapper.testMatchedAnnotations("voriconazole", PrescribingGuidanceSource.CPIC_GUIDELINE, 2);
    testWrapper.testMatchedAnnotations("voriconazole", PrescribingGuidanceSource.DPWG_GUIDELINE, 1);
  }

  @Test
  void testCyp2c19noCall(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("CYP2C19", "rs12769205", "A", "G")
        .variation("CYP2C19", "rs4244285", "A", "A");
    testWrapper.executeWithOutsideCalls(s_otherOutsideCallFilePath);

    testWrapper.testNotCalledByMatcher("CYP2C19");

    testWrapper.testNoMatchFromSource("citalopram", PrescribingGuidanceSource.CPIC_GUIDELINE);
    testWrapper.testNoMatchFromSource("citalopram", PrescribingGuidanceSource.DPWG_GUIDELINE);
    testWrapper.testNoMatchFromSource("ivacaftor", PrescribingGuidanceSource.CPIC_GUIDELINE);
    testWrapper.testNoMatchFromSource("ivacaftor", PrescribingGuidanceSource.DPWG_GUIDELINE);
  }

  @Test
  void testCyp2c19s4bs17rs28399504missing(TestInfo testInfo) throws Exception {
    // NOTE: this test has multiple annotations for a single population - amitriptyline
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("CYP2C19", "rs12248560", "T", "T")
        .missing("CYP2C19", "rs28399504")
        .variation("CYP2C19", "rs3758581", "G", "G");
    testWrapper.execute();

    testWrapper.testCalledByMatcher("CYP2C19");
    testWrapper.testPrintCpicCalls("CYP2C19", "*4/*4", "*4/*17", "*17/*17");

    testWrapper.testMatchedAnnotations("citalopram", 8);
    testWrapper.testMatchedAnnotations("citalopram", PrescribingGuidanceSource.CPIC_GUIDELINE, 3);
    testWrapper.testMatchedAnnotations("citalopram", PrescribingGuidanceSource.DPWG_GUIDELINE, 3);
    testWrapper.testMatchedAnnotations("citalopram", PrescribingGuidanceSource.FDA_LABEL, 1);
    testWrapper.testMatchedAnnotations("citalopram", PrescribingGuidanceSource.FDA_ASSOC, 1);
  }

  @Test
  void testCyp2c19s1s4het(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("CYP2C19", "rs12248560", "T", "T")
        .variation("CYP2C19", "rs28399504", "A", "G")
        .variation("CYP2C19", "rs3758581", "G", "G");
    testWrapper.executeWithOutsideCalls(s_outsideCallFilePath);

    testWrapper.testCalledByMatcher("CYP2C19");
    testWrapper.testReportable("CYP2D6");

    testWrapper.testPrintCpicCalls("CYP2D6", "*1/*4");
    testWrapper.testPrintCpicCalls("CYP2C19", "*4/*17");

    GeneReport geneReport = testWrapper.getContext().getGeneReport(DataSource.CPIC, "CYP2D6");
    assertNotNull(geneReport);
    assertTrue(geneReport.isOutsideCall());
  }

  @Test
  void testCyp2c19s1s4missingS1(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("CYP2C19", "rs12248560", "C", "T")
        .variation("CYP2C19", "rs28399504", "A", "G")
        .missing("CYP2C19", "rs3758581");
    testWrapper.executeWithOutsideCalls(s_outsideCallFilePath);

    testWrapper.testCalledByMatcher("CYP2C19");
    testWrapper.testReportable("CYP2D6");

    testWrapper.testPrintCpicCalls("CYP2D6", "*1/*4");
    testWrapper.testPrintCpicCalls("CYP2C19",  "*1/*4", "*4/*38");

    GeneReport cyp2d6Report = testWrapper.getContext().getGeneReport(DataSource.CPIC, "CYP2D6");
    assertNotNull(cyp2d6Report);
    assertTrue(cyp2d6Report.isOutsideCall());

    GeneReport cyp2c19report = testWrapper.getContext().getGeneReport(DataSource.CPIC, "CYP2C19");
    assertNotNull(cyp2c19report);
    assertTrue(cyp2c19report.isMissingVariants());

    assertFalse(cyp2c19report.isPhased());
    assertTrue(cyp2c19report.findVariantReport("rs12248560").map(VariantReport::isHetCall).orElse(false));
    assertTrue(cyp2c19report.findVariantReport("rs3758581").map(VariantReport::isMissing).orElse(false));

    assertTrue(cyp2c19report.hasHaplotype("*38"));

    // message is for *1/*4 being ambiguous with unphased data
    assertEquals(2, cyp2c19report.getMessages().stream()
        .filter(m -> m.getExceptionType().equals(MessageAnnotation.TYPE_AMBIGUITY))
        .count());

    testWrapper.testAnyMatchFromSource("amitriptyline", PrescribingGuidanceSource.CPIC_GUIDELINE);
    testWrapper.testMessageCountForDrug(PrescribingGuidanceSource.CPIC_GUIDELINE, "amitriptyline", 2);
  }

  @Test
  void testCyp2c19SingleGeneMatch(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("CYP2C19")
        .variation("CYP2C19", "rs3758581", "A", "G")
        .missing("CYP2C19", "rs56337013");
    testWrapper.executeWithOutsideCalls(s_outsideCallFilePath);

    testWrapper.testCalledByMatcher("CYP2C19");
    testWrapper.testReportable("CYP2D6");

    testWrapper.testPrintCpicCalls("CYP2D6", "*1/*4");
    testWrapper.testPrintCpicCalls("CYP2C19", "*1/*38");

    GeneReport cyp2d6Report = testWrapper.getContext().getGeneReport(DataSource.CPIC, "CYP2D6");
    assertNotNull(cyp2d6Report);
    assertTrue(cyp2d6Report.isOutsideCall());
  }


  @Test
  void testRosuvastatin(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("ABCG2", "rs2231142", "G", "T")
        .variation("SLCO1B1", "rs56101265", "T", "C");
    Path vcfFile = testWrapper.execute();

    testWrapper.testCalledByMatcher("ABCG2", "SLCO1B1");
    testWrapper.testPrintCpicCalls("SLCO1B1", "*1/*2");

    testWrapper.testMatchedAnnotations("rosuvastatin", 1);

    // no dpyd - should not have DPYD warning
    Document document = readHtmlReport(vcfFile);
    Elements capecitabineSection = document.getElementsByClass("capecitabine");
    assertEquals(0, capecitabineSection.size());

    Elements dpydSection = document.select(".gene.dpyd");
    assertEquals(1, dpydSection.size());
    assertEquals(1, dpydSection.get(0).getElementsByClass("no-data").size());
  }

  @Test
  void testSlco1b1HomWild(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("SLCO1B1");
    Path vcfFile = testWrapper.execute();

    List<String> expectedCalls = List.of("*1/*1");

    testWrapper.testCalledByMatcher("SLCO1B1");
    testWrapper.testSourceDiplotypes(DataSource.CPIC, "SLCO1B1", expectedCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "SLCO1B1", expectedCallsToRecommendedDiplotypes(expectedCalls));
    testWrapper.testPrintCalls(DataSource.CPIC, "SLCO1B1", expectedCalls);

    testWrapper.testMatchedAnnotations("simvastatin", 2);

    DrugReport dpwgReport = testWrapper.getContext().getDrugReport(PrescribingGuidanceSource.DPWG_GUIDELINE, "simvastatin");
    assertNotNull(dpwgReport);
    assertNotNull(dpwgReport.getGuidelines());
    assertEquals("No recommendation", dpwgReport.getGuidelines().first().getAnnotations().first().getClassification());

    Document document = readHtmlReport(vcfFile);
    htmlChecks(document, "SLCO1B1", expectedCalls, "simvastatin", RecPresence.YES, RecPresence.YES);
  }

  @Test
  void testSlco1b1HomVar(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("SLCO1B1", "rs2306283", "A", "G")
        .variation("SLCO1B1", "rs4149056", "C", "C");
    Path vcfFile = testWrapper.execute();

    List<String> expectedCalls = List.of("*5/*15");

    testWrapper.testCalledByMatcher("SLCO1B1");
    testWrapper.testSourceDiplotypes(DataSource.CPIC, "SLCO1B1", expectedCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "SLCO1B1", expectedCallsToRecommendedDiplotypes(expectedCalls));
    testWrapper.testPrintCalls(DataSource.CPIC, "SLCO1B1", expectedCalls);

    testWrapper.testMatchedAnnotations("simvastatin", PrescribingGuidanceSource.CPIC_GUIDELINE, 1);
    testWrapper.testMatchedAnnotations("simvastatin", PrescribingGuidanceSource.DPWG_GUIDELINE, 1);

    Document document = readHtmlReport(vcfFile);
    htmlChecks(document, "SLCO1B1", expectedCalls, "simvastatin", RecPresence.YES, RecPresence.YES);
  }

  @Test
  void testSlco1b1Test5(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("SLCO1B1", "rs2306283", "A", "G")
        .variation("SLCO1B1", "rs11045852", "A", "G")
        .variation("SLCO1B1", "rs74064213", "A", "G");
    Path vcfFile = testWrapper.execute();

    List<String> expectedCalls = List.of("*1/*44");

    testWrapper.testCalledByMatcher("SLCO1B1");
    testWrapper.testSourceDiplotypes(DataSource.CPIC, "SLCO1B1", expectedCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "SLCO1B1", expectedCallsToRecommendedDiplotypes(expectedCalls));
    testWrapper.testPrintCalls(DataSource.CPIC, "SLCO1B1", expectedCalls);

    testWrapper.testMatchedAnnotations("simvastatin", 1);

    Document document = readHtmlReport(vcfFile);
    htmlChecks(document, "SLCO1B1", expectedCalls, "simvastatin", RecPresence.YES, RecPresence.YES_NO_MATCH);
  }

  @Test
  void testSlco1b1Test3(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("SLCO1B1", "rs2306283", "A", "G")
        .variation("SLCO1B1", "rs4149056", "T", "C");
    Path vcfFile = testWrapper.execute();

    List<String> expectedCalls = List.of("*1/*15");

    testWrapper.testCalledByMatcher("SLCO1B1");
    testWrapper.testSourceDiplotypes(DataSource.CPIC, "SLCO1B1", expectedCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "SLCO1B1", expectedCallsToRecommendedDiplotypes(expectedCalls));
    testWrapper.testPrintCalls(DataSource.CPIC, "SLCO1B1", expectedCalls);

    testWrapper.testMatchedAnnotations("simvastatin", PrescribingGuidanceSource.CPIC_GUIDELINE, 1);
    testWrapper.testMatchedAnnotations("simvastatin", PrescribingGuidanceSource.DPWG_GUIDELINE, 1);

    Document document = readHtmlReport(vcfFile);
    htmlChecks(document, "SLCO1B1", expectedCalls, "simvastatin", RecPresence.YES, RecPresence.YES);
  }

  @Test
  void testSlco1b1Test4(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("SLCO1B1", "rs4149056", "T", "C")
        .variation("SLCO1B1", "rs71581941", "C", "T");
    Path vcfFile = testWrapper.execute();

    List<String> expectedCalls = List.of("*5/*45");

    testWrapper.testCalledByMatcher("SLCO1B1");
    testWrapper.testSourceDiplotypes(DataSource.CPIC, "SLCO1B1", expectedCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "SLCO1B1", expectedCallsToRecommendedDiplotypes(expectedCalls));
    testWrapper.testPrintCalls(DataSource.CPIC, "SLCO1B1", expectedCalls);

    testWrapper.testMatchedAnnotations("simvastatin", 2);

    Document document = readHtmlReport(vcfFile);
    htmlChecks(document, "SLCO1B1", expectedCalls, "simvastatin", RecPresence.YES, RecPresence.YES_NO_MATCH);
  }

  /**
   * This tests a special case of SLCO1B1. The gene in this scenario is "uncalled" by the matcher due the sample VCF
   * data. However, SLCO1B1 has an override that will display the rs4149056 diplotype regardless of call status. That
   * same override will assign alleles to use for recommendation lookup
   */
  @Test
  void testSlco1b1UncalledOverride(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("SLCO1B1", "rs2306283", "G", "G")
        .variation("SLCO1B1", "rs4149056", "T", "C")
        .variation("SLCO1B1", "rs11045853", "A", "A")
        .variation("SLCO1B1", "rs72559748", "G", "G");
    Path vcfFile = testWrapper.execute();

    List<String> expectedCalls = List.of("rs4149056 C/rs4149056 T");

    testWrapper.testNotCalledByMatcher("SLCO1B1");
    testWrapper.testSourceDiplotypes(DataSource.CPIC, "SLCO1B1", UNKNOWN_CALL);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "SLCO1B1", List.of("*1", "*5"));
    testWrapper.testPrintCalls(DataSource.CPIC, "SLCO1B1", expectedCalls);

    testWrapper.testMatchedAnnotations("simvastatin", 3);
    testWrapper.testMatchedAnnotations("simvastatin", PrescribingGuidanceSource.CPIC_GUIDELINE, 1);
    testWrapper.testMatchedAnnotations("simvastatin", PrescribingGuidanceSource.DPWG_GUIDELINE, 1);
    testWrapper.testMatchedAnnotations("simvastatin", PrescribingGuidanceSource.FDA_ASSOC, 1);

    Document document = readHtmlReport(vcfFile);
    htmlChecks(document, "SLCO1B1", expectedCalls, "simvastatin", RecPresence.YES, RecPresence.YES);
  }


  @Test
  void testUgt1a1Phased(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .phased()
        .variation("UGT1A1", "rs887829", "C", "T");
    testWrapper.execute();

    testWrapper.testCalledByMatcher("UGT1A1");
    testWrapper.testPrintCpicCalls("UGT1A1", "*1/*80");
    testWrapper.testRecommendedDiplotypes("UGT1A1", "*1", "*80");
  }

  @Test
  void testUgt1a1Unphased(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .phased()
        .variation("UGT1A1", "rs887829", "C", "T");
    testWrapper.execute();

    testWrapper.testCalledByMatcher("UGT1A1");
    testWrapper.testPrintCpicCalls("UGT1A1", "*1/*80");
    testWrapper.testRecommendedDiplotypes("UGT1A1", "*1", "*80");
  }

  @Test
  void testUgt1a1s1s1(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("UGT1A1");
    testWrapper.execute();

    testWrapper.testCalledByMatcher("UGT1A1");
    testWrapper.testPrintCpicCalls("UGT1A1", "*1/*1");
    testWrapper.testRecommendedDiplotypes("UGT1A1", "*1", "*1");
  }

  @Test
  void testUgt1a1S1S80S28(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("UGT1A1", "rs887829", "C", "T")
        .variation("UGT1A1", "rs3064744", "TA(7)", "TA(8)");
    testWrapper.execute();

    testWrapper.testCalledByMatcher("UGT1A1");
    testWrapper.testPrintCpicCalls("UGT1A1", "*1/*80+*28");
    testWrapper.testRecommendedDiplotypes("UGT1A1", "*1", "*80+*28");
  }

  @Test
  void testUgt1a1S28S37(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("UGT1A1", "rs3064744", "TA(8)", "TA(9)");
    testWrapper.execute();

    testWrapper.testCalledByMatcher("UGT1A1");
    testWrapper.testPrintCpicCalls("UGT1A1", "*28/*37");
    testWrapper.testRecommendedDiplotypes("UGT1A1", "*37", "*28");
  }

  @Test
  void testUgt1a1s28s80phased(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .phased()
        .variation("UGT1A1", "rs887829", "C", "T")
        .variation("UGT1A1", "rs3064744", "TA(7)", "TA(8)");
    testWrapper.execute();

    testWrapper.testCalledByMatcher("UGT1A1");
    testWrapper.testPrintCpicCalls("UGT1A1", "*1/*80+*28");
    testWrapper.testRecommendedDiplotypes("UGT1A1", "*1", "*80+*28");

    // the guideline should not have an ambiguity message
    testWrapper.testMatchedAnnotations("atazanavir", 1);
    testWrapper.testMessageCountForDrug(PrescribingGuidanceSource.CPIC_GUIDELINE, "atazanavir", 0);

    testWrapper.testMessageCountForGene(DataSource.CPIC, "UGT1A1", 2);
    testWrapper.testGeneHasMessage(DataSource.CPIC, "UGT1A1", "reference-allele");
  }

  @Test
  void testUgt1a1s28s80s6s60phased(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .phased()
        .variation("UGT1A1", "rs887829", "T", "C")
        .variation("UGT1A1", "rs3064744", "TA(8)", "TA(7)")
        .variation("UGT1A1", "rs4148323", "G", "A");
    testWrapper.execute();

    testWrapper.testCalledByMatcher("UGT1A1");
    testWrapper.testPrintCpicCalls("UGT1A1", "*6/*80+*28");
    testWrapper.testRecommendedDiplotypes("UGT1A1", "*6", "*80+*28");
  }

  @Test
  void testUgt1a1s28s80s6s60unphased(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("UGT1A1", "rs887829", "T", "C")
        .variation("UGT1A1", "rs3064744", "TA(8)", "TA(7)")
        .variation("UGT1A1", "rs4148323", "G", "A");
    testWrapper.execute();

    testWrapper.testCalledByMatcher("UGT1A1");
    testWrapper.testPrintCpicCalls("UGT1A1", "*6/*80+*28");
    testWrapper.testRecommendedDiplotypes("UGT1A1", "*6", "*80+*28");
  }

  @Test
  void testUgt1a1s6s6(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("UGT1A1", "rs4148323", "A", "A");
    testWrapper.execute();

    testWrapper.testCalledByMatcher("UGT1A1");
    testWrapper.testPrintCpicCalls("UGT1A1", "*6/*6");
    testWrapper.testRecommendedDiplotypes("UGT1A1", "*6", "*6");
  }

  @Test
  void testUgt1a1s6s60s80s28MissingPhased(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .phased()
        .missing("UGT1A1", "rs3064744")
        .variation("UGT1A1", "rs887829", "C", "T")
        .variation("UGT1A1", "rs4148323", "A", "G");
    testWrapper.execute();

    testWrapper.testCalledByMatcher("UGT1A1");
    testWrapper.testPrintCpicCalls("UGT1A1", "*6/*80", "*6/*80+*28", "*6/*80+*37");
    testWrapper.testRecommendedDiplotypes("UGT1A1", "*6", "*80");
    testWrapper.testRecommendedDiplotypes("UGT1A1", "*6", "*80+*28");
    testWrapper.testRecommendedDiplotypes("UGT1A1", "*6", "*80+*37");
  }

  @Test
  void testUgt1a1s6s60s80s28MissingUnphased(TestInfo testInfo) throws Exception {
    // NOTE: this test has multiple annotations for a single population - atazanavir
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .missing("UGT1A1", "rs3064744")
        .variation("UGT1A1", "rs887829", "C", "T")
        .variation("UGT1A1", "rs4148323", "A", "G");
    testWrapper.execute();

    testWrapper.testCalledByMatcher("UGT1A1");
    testWrapper.testPrintCpicCalls("UGT1A1", "*6/*80", "*6/*80+*28", "*6/*80+*37");
    testWrapper.testRecommendedDiplotypes("UGT1A1", "*6", "*80");
    testWrapper.testRecommendedDiplotypes("UGT1A1", "*6", "*80+*28");
    testWrapper.testRecommendedDiplotypes("UGT1A1", "*6", "*80+*37");
  }

  @Test
  void testUgt1a1s80s28missing(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .missing("UGT1A1", "rs3064744")
        .variation("UGT1A1", "rs887829", "C", "T");
    Path vcfFile = testWrapper.execute();

    testWrapper.testCalledByMatcher("UGT1A1");
    testWrapper.testRecommendedDiplotypes("UGT1A1", "*1", "*80");
    testWrapper.testRecommendedDiplotypes("UGT1A1", "*1", "*80+*28");
    testWrapper.testRecommendedDiplotypes("UGT1A1", "*1", "*80+*37");
    testWrapper.testPrintCpicCalls("UGT1A1", "*1/*80", "*1/*80+*28", "*1/*80+*37");

    Document document = readHtmlReport(vcfFile);
    Elements drugSection = document.select(".cpic-guideline-atazanavir");
    assertEquals(2, drugSection.size());

    Elements d1 = drugSection.get(0).select(".rx-dip");
    assertEquals(1, d1.size());
    assertEquals("UGT1A1:*1/*80", d1.get(0).text());

    Elements d2 = drugSection.get(1).select(".rx-dip");
    assertEquals(2, d2.size());
    assertEquals(List.of("UGT1A1:*1/*80+*28", "UGT1A1:*1/*80+*37"), d2.stream().map(Element::text).toList());
  }

  @Test
  void testUgt1a1na12717(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("UGT1A1", "rs887829", "T", "T")
        .variation("UGT1A1", "rs3064744", "TA(8)", "TA(7)");
    testWrapper.execute();

    testWrapper.testCalledByMatcher("UGT1A1");
    testWrapper.testPrintCpicCalls("UGT1A1", "*80/*80+*28");
    testWrapper.testRecommendedDiplotypes("UGT1A1", "*80", "*80+*28");
  }

  @Test
  void testUgt1a1s28homMissing(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .missing("UGT1A1", "rs887829")
        .variation("UGT1A1", "rs3064744", "TA(8)", "TA(8)");
    testWrapper.execute();

    testWrapper.testCalledByMatcher("UGT1A1");
    testWrapper.testPrintCpicCalls("UGT1A1", "*28/*28", "*28/*80+*28", "*80+*28/*80+*28");
    testWrapper.testRecommendedDiplotypes("UGT1A1", "*28", "*28");
    testWrapper.testRecommendedDiplotypes("UGT1A1", "*28", "*80+*28");
    testWrapper.testRecommendedDiplotypes("UGT1A1", "*80+*28", "*80+*28");
  }

  @Test
  void testUgt1a1s28s60Hom(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("UGT1A1", "rs3064744", "TA(8)", "TA(7)");
    testWrapper.execute();

    testWrapper.testCalledByMatcher("UGT1A1");
    testWrapper.testPrintCpicCalls("UGT1A1", "*1/*28");
    testWrapper.testRecommendedDiplotypes("UGT1A1", "*1", "*28");
  }

  @Test
  void testUgt1a1s27s28unphaseds80s60missing(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .missing("UGT1A1", "rs887829")
        .variation("UGT1A1", "rs3064744", "TA(8)", "TA(7)")
        .variation("UGT1A1", "rs35350960", "C", "A");
    testWrapper.execute();

    testWrapper.testCalledByMatcher("UGT1A1");
    testWrapper.testPrintCpicCalls("UGT1A1", "*27/*28", "*27/*80+*28");
    testWrapper.testRecommendedDiplotypes("UGT1A1", "*27", "*28");
    testWrapper.testRecommendedDiplotypes("UGT1A1", "*27", "*80+*28");
  }

  @Test
  void testUgt1a1HG00436(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .phased()
        .variation("UGT1A1", "rs887829", "T", "C")
        .variation("UGT1A1", "rs3064744", "TA(8)", "TA(7)")
        .variation("UGT1A1", "rs35350960", "A", "C");
    testWrapper.execute();

    testWrapper.testNotCalledByMatcher("UGT1A1");
  }

  @Test
  void testUgt1a1s1s80s27s60s28MissingPhased(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .phased()
        .missing("UGT1A1", "rs3064744")
        .variation("UGT1A1", "rs887829", "T", "C")
        .variation("UGT1A1", "rs35350960", "A", "C");
    testWrapper.execute();

    testWrapper.testNotCalledByMatcher("UGT1A1");
    GeneReport geneReport = testWrapper.getContext().getGeneReport(DataSource.CPIC, "UGT1A1");
    assertNotNull(geneReport);
    assertTrue(geneReport.isPhased());
  }

  @Test
  void testUgt1a1s1s60s80s6phased(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .phased()
        .variation("UGT1A1", "rs887829", "T", "C")
        .variation("UGT1A1", "rs35350960", "A", "C");
    testWrapper.execute();

    testWrapper.testNotCalledByMatcher("UGT1A1");
    GeneReport geneReport = testWrapper.getContext().getGeneReport(DataSource.CPIC, "UGT1A1");
    assertNotNull(geneReport);
    assertTrue(geneReport.isPhased());
  }

  @Test
  void testUgt1a1s1s60s80s28s6phased(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .phased()
        .variation("UGT1A1", "rs887829", "T", "C")
        .variation("UGT1A1", "rs3064744", "TA(8)", "TA(7)")
        .variation("UGT1A1", "rs35350960", "A", "C");
    testWrapper.execute();

    testWrapper.testNotCalledByMatcher("UGT1A1");
    GeneReport geneReport = testWrapper.getContext().getGeneReport(DataSource.CPIC, "UGT1A1");
    assertNotNull(geneReport);
    assertTrue(geneReport.isPhased());
  }

  @Test
  void testUgt1a1s1s37s80s60phased(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .phased()
        .variation("UGT1A1", "rs887829", "T", "C")
        .variation("UGT1A1", "rs3064744", "TA(9)", "TA(7)");
    testWrapper.execute();

    testWrapper.testCalledByMatcher("UGT1A1");
    testWrapper.testReportable("UGT1A1");
    testWrapper.testPrintCpicCalls("UGT1A1", "*1/*80+*37");
    testWrapper.testRecommendedDiplotypes("UGT1A1", "*1", "*80+*37");
    GeneReport geneReport = testWrapper.getContext().getGeneReport(DataSource.CPIC, "UGT1A1");
    assertNotNull(geneReport);
    assertTrue(geneReport.isPhased());
  }

  @Test
  void testCyp3a5Missing3Message(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .missing("CYP3A5", "rs776746");
    testWrapper.execute();

    testWrapper.testCalledByMatcher("CYP3A5");
    testWrapper.testReportable("CYP3A5");
    testWrapper.testPrintCpicCalls("CYP3A5", "*1/*1");
    testWrapper.testRecommendedDiplotypes("CYP3A5", "*1", "*1");

    GeneReport gene = testWrapper.getContext().getGeneReport(DataSource.CPIC, "CYP3A5");
    assertNotNull(gene);
    // rs776746 should be missing from this report
    assertNotNull(gene.getVariantReports());
    assertTrue(gene.getVariantReports().stream().anyMatch(v -> v.isMissing() && v.getDbSnpId().equals("rs776746")));

    // the guideline should have a matching message
    assertTrue(testWrapper.getContext().getDrugReports().get(PrescribingGuidanceSource.CPIC_GUIDELINE).values().stream()
        .filter(r -> r.getName().equals("tacrolimus"))
        .noneMatch(r -> r.getMessages().isEmpty()));

    assertFalse(gene.isPhased());
  }

  @Test
  void testCyp3a5v1(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("CYP3A5", "rs776746", "T", "C");
    testWrapper.execute();

    testWrapper.testCalledByMatcher("CYP3A5");
    testWrapper.testReportable("CYP3A5");
    testWrapper.testPrintCpicCalls("CYP3A5", "*1/*3");
    testWrapper.testRecommendedDiplotypes("CYP3A5", "*1", "*3");
  }

  @Test
  void testCyp3a5v2(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("CYP3A5", "rs28383479", "C", "T")
        .variation("CYP3A5", "rs776746", "C", "T")
    ;
    testWrapper.execute();

    testWrapper.testCalledByMatcher("CYP3A5");
    testWrapper.testPrintCpicCalls("CYP3A5", "*3/*9");
    testWrapper.testRecommendedDiplotypes("CYP3A5", "*3", "*9");
  }

  @Test
  void testCyp3a5v3(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("CYP3A5", "rs776746", "C", "C")
    ;
    testWrapper.execute();

    testWrapper.testCalledByMatcher("CYP3A5");
    testWrapper.testPrintCpicCalls("CYP3A5", "*3/*3");
    testWrapper.testRecommendedDiplotypes("CYP3A5", "*3", "*3");
  }

  @Test
  void testCyp3a5v4(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("CYP3A5", "rs776746", "T", "C")
    ;
    testWrapper.execute();

    testWrapper.testCalledByMatcher("CYP3A5");
    testWrapper.testPrintCpicCalls("CYP3A5", "*1/*3");
    testWrapper.testRecommendedDiplotypes("CYP3A5", "*1", "*3");
  }

  @Test
  void testCyp3a5v5(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("CYP3A5", "rs28383479", "T", "C")
        .variation("CYP3A5", "rs776746", "T", "C")
    ;
    testWrapper.execute();

    testWrapper.testCalledByMatcher("CYP3A5");
    testWrapper.testPrintCpicCalls("CYP3A5", "*3/*9");
    testWrapper.testRecommendedDiplotypes("CYP3A5", "*3", "*9");
  }

  @Test
  void testHlab(TestInfo testInfo) throws Exception {
    Path outsideCallPath = TestUtils.createTestFile(testInfo, ".tsv");
    try (PrintWriter writer = new PrintWriter(Files.newBufferedWriter(outsideCallPath))) {
      writer.println("HLA-B\t*15:02/*57:01");
    }

    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("CYP2C9");
    testWrapper.executeWithOutsideCalls(outsideCallPath);

    testWrapper.testCalledByMatcher("CYP2C9");
    testWrapper.testReportable("CYP2C9");
    testWrapper.testPrintCalls(DataSource.CPIC, "CYP2C9", "*1/*1");
    testWrapper.testPrintCalls(DataSource.DPWG, "CYP2C9", "*1/*1");

    testWrapper.testNotCalledByMatcher("HLA-B");
    testWrapper.testReportable("HLA-B");
    testWrapper.testSourcePhenotype(DataSource.CPIC, "HLA-B", "*57:01 positive");
    testWrapper.testSourcePhenotype(DataSource.CPIC, "HLA-B", "*58:01 negative");
    testWrapper.testSourcePhenotype(DataSource.CPIC, "HLA-B", "*15:02 positive");
    testWrapper.testSourcePhenotype(DataSource.DPWG, "HLA-B", "*57:01 positive");
    testWrapper.testSourcePhenotype(DataSource.DPWG, "HLA-B", "*58:01 negative");
    testWrapper.testSourcePhenotype(DataSource.DPWG, "HLA-B", "*15:02 positive");

    // *57:01 guideline
    testWrapper.testMatchedAnnotations("abacavir", PrescribingGuidanceSource.CPIC_GUIDELINE, 1);
    testWrapper.testMatchedAnnotations("abacavir", PrescribingGuidanceSource.DPWG_GUIDELINE, 1);
    // *58:01 guideline
    testWrapper.testMatchedAnnotations("allopurinol", PrescribingGuidanceSource.CPIC_GUIDELINE, 1);
    // *15:02 guideline (along with CYP2C9)
    testWrapper.testMatchedAnnotations("phenytoin", 6);
    testWrapper.testMatchedAnnotations("phenytoin", PrescribingGuidanceSource.CPIC_GUIDELINE, 2);
    testWrapper.testMatchedAnnotations("phenytoin", PrescribingGuidanceSource.DPWG_GUIDELINE, 2);
    testWrapper.testMatchedAnnotations("phenytoin", PrescribingGuidanceSource.FDA_LABEL, 1);
    testWrapper.testMatchedAnnotations("phenytoin", PrescribingGuidanceSource.FDA_ASSOC, 1);
  }

  @Test
  void testSingleHlabAllele(TestInfo testInfo) throws Exception {
    Path outsideCallPath = TestUtils.createTestFile(testInfo, ".tsv");
    try (PrintWriter writer = new PrintWriter(Files.newBufferedWriter(outsideCallPath))) {
      writer.println("HLA-B\t*15:02");
    }

    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("CYP2C9");
    testWrapper.executeWithOutsideCalls(outsideCallPath);

    testWrapper.testNotCalledByMatcher("HLA-B");
    testWrapper.testReportable("HLA-B");
    testWrapper.testSourcePhenotype(DataSource.CPIC, "HLA-B", "*57:01 negative");
    testWrapper.testSourcePhenotype(DataSource.CPIC, "HLA-B", "*58:01 negative");
    testWrapper.testSourcePhenotype(DataSource.CPIC, "HLA-B", "*15:02 positive");
    testWrapper.testSourcePhenotype(DataSource.DPWG, "HLA-B", "*57:01 negative");
    testWrapper.testSourcePhenotype(DataSource.DPWG, "HLA-B", "*58:01 negative");
    testWrapper.testSourcePhenotype(DataSource.DPWG, "HLA-B", "*15:02 positive");

    testWrapper.testMatchedAnnotations("abacavir", PrescribingGuidanceSource.CPIC_GUIDELINE, 1);
    testWrapper.testMatchedAnnotations("allopurinol", PrescribingGuidanceSource.CPIC_GUIDELINE, 1);
    testWrapper.testMatchedAnnotations("phenytoin", 6);
    testWrapper.testMatchedAnnotations("phenytoin", PrescribingGuidanceSource.CPIC_GUIDELINE, 2);
    testWrapper.testMatchedAnnotations("phenytoin", PrescribingGuidanceSource.DPWG_GUIDELINE, 2);
    testWrapper.testMatchedAnnotations("phenytoin", PrescribingGuidanceSource.FDA_LABEL, 1);
    testWrapper.testMatchedAnnotations("phenytoin", PrescribingGuidanceSource.FDA_ASSOC, 1);

    // carbamazepine-CPIC is a two gene lookup, let's test to make sure only specifying HLA-B will still return results
    testWrapper.testMatchedAnnotations("carbamazepine", PrescribingGuidanceSource.CPIC_GUIDELINE, 3);
    // carbamazepine-DPWG is a single gene lookup so it will be a "normal" lookup
    testWrapper.testMatchedAnnotations("carbamazepine", PrescribingGuidanceSource.DPWG_GUIDELINE, 1);
  }

  @Test
  void testHlabPhenotype(TestInfo testInfo) throws Exception {
    Path outsideCallPath = TestUtils.createTestFile(testInfo, ".tsv");
    try (PrintWriter writer = new PrintWriter(Files.newBufferedWriter(outsideCallPath))) {
      writer.println("HLA-B\t\t*57:01 positive");
    }

    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("CYP2C9");
    testWrapper.executeWithOutsideCalls(outsideCallPath);

    testWrapper.testCalledByMatcher("CYP2C9");
    testWrapper.testNotCalledByMatcher("HLA-B");
    testWrapper.testReportable("CYP2C9");
    testWrapper.testReportable("HLA-B");
    testWrapper.testMatchedAnnotations("abacavir", 4);
    testWrapper.testMatchedAnnotations("abacavir", PrescribingGuidanceSource.CPIC_GUIDELINE, 1);
    testWrapper.testMatchedAnnotations("abacavir", PrescribingGuidanceSource.DPWG_GUIDELINE, 1);
    testWrapper.testMatchedAnnotations("abacavir", PrescribingGuidanceSource.FDA_LABEL, 1);
    testWrapper.testMatchedAnnotations("abacavir", PrescribingGuidanceSource.FDA_ASSOC, 1);
    // allopurinol relies on a different allele for recs so no matches
    testWrapper.testMatchedAnnotations("allopurinol", 0);
    // phenytoin also relies on a different allele, but there will be a match for DPWG since the recommendations are
    // split between the two genes on that side
    testWrapper.testMatchedAnnotations("phenytoin", 1);
    testWrapper.testNoMatchFromSource("phenytoin", PrescribingGuidanceSource.CPIC_GUIDELINE);
    testWrapper.testAnyMatchFromSource("phenytoin", PrescribingGuidanceSource.DPWG_GUIDELINE);
  }

  /**
   * An example report that shows a few different types of recommendation scenarios all in one report. The examples
   * shown are:
   * <ul>
   *   <li>celecoxib = 1 CPIC recommendation</li>
   *   <li>citalopram = 2 recommendations: 1 CPIC, 1 DPWG, 1 gene and it's called</li>
   *   <li>clomipramine = 2 recommendations: 1 CPIC, 1 DPWG, 2 gene but only 1 called</li>
   *   <li>carbamazepine = 3 CPIC recommendations on different populations</li>
   *   <li>clopidogrel = 4 recommendations: 3 CPIC on different pops, 1 DPWG</li>
   *   <li>flucloxacillin = 1 recommendation from DPWG but none exists for CPIC</li>
   *   <li>fluvoxamine = 0 recommendations, no gene reportable</li>
   *   <li>siponimod = 1 DPWG recommendation</li>
   * </ul>
   */
  @Test
  void testRecommendationExamples(TestInfo testInfo) throws Exception {
    Path outsideCallPath = TestUtils.createTestFile(testInfo, ".tsv");
    try (PrintWriter writer = new PrintWriter(Files.newBufferedWriter(outsideCallPath))) {
      writer.println("HLA-A\t\t*31:01 positive");
      writer.println("HLA-B\t*57:01/*58:01\t");
    }

    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("CYP2C9")
        .variation("CYP2C19", "rs12769205", "G", "G")
        .variation("CYP2C19", "rs4244285", "A", "A")
        .variation("CYP2C19", "rs3758581", "G", "G");
    testWrapper.executeWithOutsideCalls(outsideCallPath);

    testWrapper.testRecommendedDiplotypes("CYP2C19", "*2", "*2");
    testWrapper.testPrintCpicCalls("CYP2C19", "*2/*2");
    testWrapper.testNotCalledByMatcher("CYP2D6");

    GeneReport cyp2c9 = testWrapper.getContext().getGeneReport(DataSource.CPIC, "CYP2C9");
    assertNotNull(cyp2c9);
    assertEquals(1, cyp2c9.getRecommendationDiplotypes().size());
    assertTrue(cyp2c9.getRecommendationDiplotypes().stream().allMatch(d -> d.getActivityScore().equals("2.0")));

    testWrapper.testReportable("CYP2C19", "CYP2C9", "HLA-A", "HLA-B");
    testWrapper.testMatchedAnnotations("celecoxib", 1);
    testWrapper.testAnyMatchFromSource("celecoxib", PrescribingGuidanceSource.CPIC_GUIDELINE);
    testWrapper.testMatchedAnnotations("citalopram", 4);
    testWrapper.testMatchedAnnotations("clomipramine", PrescribingGuidanceSource.CPIC_GUIDELINE, 1);
    testWrapper.testMatchedAnnotations("clomipramine", PrescribingGuidanceSource.DPWG_GUIDELINE, 1);
    testWrapper.testMatchedAnnotations("clopidogrel", 6);
    testWrapper.testMatchedAnnotations("clopidogrel", PrescribingGuidanceSource.CPIC_GUIDELINE, 3);
    testWrapper.testMatchedAnnotations("clopidogrel", PrescribingGuidanceSource.DPWG_GUIDELINE, 1);
    testWrapper.testMatchedAnnotations("clopidogrel", PrescribingGuidanceSource.FDA_LABEL, 1);
    testWrapper.testMatchedAnnotations("clopidogrel", PrescribingGuidanceSource.FDA_ASSOC, 1);
    testWrapper.testNoMatchFromSource("flucloxacillin", PrescribingGuidanceSource.CPIC_GUIDELINE);
    testWrapper.testMatchedAnnotations("flucloxacillin", PrescribingGuidanceSource.DPWG_GUIDELINE, 1);
    testWrapper.testNoMatchFromSource("fluvoxamine", PrescribingGuidanceSource.CPIC_GUIDELINE);
    testWrapper.testNoMatchFromSource("fluvoxamine", PrescribingGuidanceSource.DPWG_GUIDELINE);

    // siponimod has DPWG & FDA recs, DPWG uses traditional matching and FDA uses diplotype-specific matching
    testWrapper.testMatchedAnnotations("siponimod", 2);
    testWrapper.testAnyMatchFromSource("siponimod", PrescribingGuidanceSource.DPWG_GUIDELINE);
    testWrapper.testAnyMatchFromSource("siponimod", PrescribingGuidanceSource.FDA_LABEL);

    testWrapper.testMatchedAnnotations("carbamazepine", PrescribingGuidanceSource.CPIC_GUIDELINE, 3);
    testWrapper.testMatchedAnnotations("carbamazepine", PrescribingGuidanceSource.DPWG_GUIDELINE, 1);
  }

  @Test
  void testTpmtStar1s(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("TPMT", "rs1800460", "C", "T")
        .variation("TPMT", "rs1142345", "T", "C");
    testWrapper.execute();

    testWrapper.testCalledByMatcher("TPMT");
    testWrapper.testPrintCpicCalls("TPMT", "*1/*3A");
    testWrapper.testRecommendedDiplotypes("TPMT", "*1", "*3A");

    GeneReport tpmtReport = testWrapper.getContext().getGeneReport(DataSource.CPIC, "TPMT");
    assertNotNull(tpmtReport);
    assertEquals(43, tpmtReport.getVariantReports().size());
  }


  @Test
  void testCyp2c9star61(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("CYP2C9", "rs1799853", "C", "T")
        .variation("CYP2C9", "rs202201137", "A", "G");
    testWrapper.execute();

    testWrapper.testCalledByMatcher("CYP2C9");
    testWrapper.testPrintCpicCalls("CYP2C9", "*1/*61");
    testWrapper.testRecommendedDiplotypes("CYP2C9", "*1", "*61");
  }

  @Test
  void testCyp2c9star1Hom(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("CYP2C9");
    testWrapper.execute();

    testWrapper.testCalledByMatcher("CYP2C9");
    testWrapper.testPrintCpicCalls("CYP2C9", "*1/*1");
    testWrapper.testRecommendedDiplotypes("CYP2C9", "*1", "*1");
    testWrapper.testMatchedAnnotations("celecoxib", 1);
    testWrapper.testMatchedAnnotations("ibuprofen", 1);
    testWrapper.testMatchedAnnotations("lornoxicam", 1);
  }


  /**
   * Test CYP2B6 for a het *34 sample file. When doing the "top match" scenario, this will only match to a *1/*34 and
   * thus only match to a single recommendation.
   * <p>
   * This test will have a different outcome when run in "all matches" mode and should be compared with
   * {@link #testCyp2b6star1star34AllMatch(TestInfo)}.
   */
  @Test
  void testCyp2b6star1star34(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("CYP2B6", "rs34223104", "T", "C")
        .variation("CYP2B6", "rs3211371", "C", "A")
        .variation("CYP2B6", "rs3745274", "G", "T")
        .variation("CYP2B6", "rs2279343", "A", "G")
    ;
    testWrapper.execute();

    testWrapper.testCalledByMatcher("CYP2B6");
    testWrapper.testPrintCpicCalls("CYP2B6", "*1/*34");
    testWrapper.testRecommendedDiplotypes("CYP2B6", "*1", "*34");
    testWrapper.testMatchedAnnotations("efavirenz", 1);
  }

  /**
   * This test is just like {@link #testCyp2b6star1star34(TestInfo)} but run in "all matches" mode. This should result
   * in 2 possible different calls coming from the matcher. These two have different phenotypes and, thus, match to
   * different recommendations.
   */
  @Test
  void testCyp2b6star1star34AllMatch(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, true);
    testWrapper.getVcfBuilder()
        .variation("CYP2B6", "rs34223104", "T", "C")
        .variation("CYP2B6", "rs3211371", "C", "A")
        .variation("CYP2B6", "rs3745274", "G", "T")
        .variation("CYP2B6", "rs2279343", "A", "G")
    ;
    testWrapper.execute();

    testWrapper.testCalledByMatcher("CYP2B6");
    testWrapper.testPrintCpicCalls("CYP2B6", "*1/*34", "*33/*36");
    testWrapper.testRecommendedDiplotypes("CYP2B6", "*1", "*34");
    testWrapper.testRecommendedDiplotypes("CYP2B6", "*33", "*36");
    testWrapper.testMatchedAnnotations("efavirenz", 2);
  }


  /* MT-RNR1 */
  @Test
  void testMtrnr1(TestInfo testInfo) throws Exception {
    Path outsideCallPath = TestUtils.createTestFile(testInfo, ".tsv");
    try (PrintWriter writer = new PrintWriter(Files.newBufferedWriter(outsideCallPath))) {
      writer.println("MT-RNR1\tm.1555A>G");
    }
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("CYP2C19")
        .reference("CYP2C9")
    ;
    testWrapper.executeWithOutsideCalls(outsideCallPath);

    testWrapper.testCalledByMatcher("CYP2C19");
    testWrapper.testCalledByMatcher("CYP2C9");
    testWrapper.testReportable("MT-RNR1");
    testWrapper.testMatchedAnnotations("amikacin", PrescribingGuidanceSource.CPIC_GUIDELINE, 1);
    testWrapper.testMatchedAnnotations("amikacin", PrescribingGuidanceSource.FDA_LABEL, 1);
  }


  /* IFNL3/4 */
  @Test
  void testIfnl3(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("IFNL3")
    ;
    testWrapper.execute();

    testWrapper.testCalledByMatcher("IFNL3");
    testWrapper.testReportable("IFNL3");
    testWrapper.testPrintCpicCalls("IFNL3", "rs12979860 reference (C)/rs12979860 reference (C)");
    testWrapper.testMatchedAnnotations("peginterferon alfa-2a", 0);
    testWrapper.testMatchedAnnotations("peginterferon alfa-2b", 0);
  }


  @Test
  void testCyp3a4(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("CYP3A4", "rs72552799", "T", "T")
    ;
    testWrapper.execute();
    testWrapper.testCalledByMatcher("CYP3A4");
    testWrapper.testReportable("CYP3A4");
    testWrapper.testPrintCalls(DataSource.DPWG, "CYP3A4", "*8/*8");
    testWrapper.testMatchedAnnotations("quetiapine", 1);
  }

  /**
   * Added to check the output of a partial match for CYP2C19 and make sure messages are applied
   */
  @Test
  void testPartialCall(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, true, true, false);
    testWrapper.getVcfBuilder()
        .reference("CYP2C19")
        .variation("CYP2C19", "rs367543002", "C", "T")
        .variation("CYP2C19", "rs3758581", "G", "G")
        .missing("CYP2C19", "rs367543003");
    testWrapper.execute();
    testWrapper.testCalledByMatcher("CYP2C19");
    testWrapper.testReportable("CYP2C19");
  }

  /**
   * Added to check the output of a partial match for CYP2C19 with a call for 2D6 to see how a two-gene recommendation
   * works with one gene being a partial call.
   */
  @Test
  void testPartialCallInTwoGene(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, true, true, false);
    testWrapper.getVcfBuilder()
        .reference("CYP2C19")
        .variation("CYP2C19", "rs367543002", "C", "T")
        .variation("CYP2C19", "rs3758581", "G", "G")
        .missing("CYP2C19", "rs367543003");
    testWrapper.executeWithOutsideCalls(s_outsideCallFilePath); //CYP2D6 *1/*4

    testWrapper.testCalledByMatcher("CYP2C19");
    testWrapper.testReportable("CYP2C19");
    testWrapper.testReportable("CYP2D6");

    // We don't expect to get a match since 2C19 is a partial call even when 2D6 has a call.
    // Currently, partials will not be visible past the phenotyper so recommendation matches are not visible.
    testWrapper.testMatchedAnnotations("amitriptyline", PrescribingGuidanceSource.CPIC_GUIDELINE, 0);
  }


  /**
   * Tests how PharmCAT handles that state when sample VCF data exists for a gene and an outside call also exists for
   * that gene. Currently, this should execute successfully by ignoring VCF data and using the outside call
   */
  @Test
  void testOutsideCallCollision(TestInfo testInfo) throws Exception {
    Path outsideCallPath = TestUtils.createTestFile(testInfo, ".tsv");
    try (PrintWriter writer = new PrintWriter(Files.newBufferedWriter(outsideCallPath))) {
      writer.println("CYP2C19\t*2/*2");
    }

    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("CYP2C19");
    testWrapper.executeWithOutsideCalls(outsideCallPath);

    testWrapper.testNotCalledByMatcher("CYP2C19");
    // this is the diplotype indicated in the outside call, not the one matched
    testWrapper.testPrintCpicCalls( "CYP2C19", "*2/*2");

    testWrapper.testMessageCountForGene(DataSource.CPIC, "CYP2C19", 2);
    testWrapper.testGeneHasMessage(DataSource.CPIC, "CYP2C19", "prefer-sample-data",
        MessageHelper.MSG_OUTSIDE_CALL);
  }



  @Test
  void outsideCallCollision2Files(TestInfo testInfo) throws Exception {

    Path outsideCallPath1 = TestUtils.createTestFile(testInfo, "1.tsv");
    try (PrintWriter writer = new PrintWriter(Files.newBufferedWriter(outsideCallPath1))) {
      writer.println("CYP4F2\t*1/*3");
    }
    Path outsideCallPath2 = TestUtils.createTestFile(testInfo, "2.tsv");
    try (PrintWriter writer = new PrintWriter(Files.newBufferedWriter(outsideCallPath2))) {
      writer.println("CYP4F2\t*9/*4");
    }

    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("CYP2C9");
    Path vcfFile = testWrapper.executeWithOutsideCalls(outsideCallPath1, outsideCallPath2);

    // this is an outside calls
    testWrapper.testNotCalledByMatcher("CYP4F2");
    // this is a regular call
    testWrapper.testCalledByMatcher("CYP2C9");

    testWrapper.testPrintCpicCalls( "CYP4F2", "*1/*3", "*4/*9");
    testWrapper.testGeneHasMessage(DataSource.CPIC, "CYP4F2", MessageHelper.MSG_OUTSIDE_CALL);

    Document document = readHtmlReport(vcfFile);
    assertEquals(1,
        document.select(".gene.CYP4F2 .alert-warning." + MessageHelper.MSG_OUTSIDE_CALL).size());
  }


  /**
   * Tests that an "unordered" diplotype should normalize to the ordered version then it can be used for matching
   */
  @Test
  void testOutsideCallDiplotypeNormalization(TestInfo testInfo) throws Exception {
    Path outsideCallPath = TestUtils.createTestFile(testInfo, ".tsv");
    try (PrintWriter writer = new PrintWriter(Files.newBufferedWriter(outsideCallPath))) {
      // diplotype in backwards order
      writer.println("CYP2C19\t*2/*1");
    }

    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("CYP2C9");
    testWrapper.executeWithOutsideCalls(outsideCallPath);

    testWrapper.testNotCalledByMatcher("CYP2C19");
    // this should be a normalized version of the given diplotype
    testWrapper.testPrintCpicCalls( "CYP2C19", "*1/*2");
  }

  @Test
  void testOutsideCallPhenotypeOverridesDiplotype(TestInfo testInfo) throws Exception {
    Path outsideCallPath = TestUtils.createTestFile(testInfo,".tsv");
    try (PrintWriter writer = new PrintWriter(Files.newBufferedWriter(outsideCallPath))) {
      writer.println("CYP2D6\t*1/*1\tPM\t" + TextConstants.GTE + "4.0");
    }

    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("CYP2C19");
    testWrapper.executeWithOutsideCalls(outsideCallPath);

    GeneReport geneReport = testWrapper.getContext().getGeneReport(DataSource.CPIC, "CYP2D6");
    assertNotNull(geneReport);
    assertEquals(1, geneReport.getRecommendationDiplotypes().size());

    Diplotype diplotype = geneReport.getRecommendationDiplotypes().first();
    assertThat(diplotype.getPhenotypes(), contains("Poor Metabolizer"));

    DrugReport drugReport = testWrapper.getContext().getDrugReport(PrescribingGuidanceSource.CPIC_GUIDELINE, "clomipramine");
    assertNotNull(drugReport);
    assertEquals(1, drugReport.getGuidelines().size());
    GuidelineReport guidelineReport = drugReport.getGuidelines().first();
    assertEquals(1, guidelineReport.getAnnotations().size());
    AnnotationReport annotationReport = guidelineReport.getAnnotations().first();
    assertEquals("Poor Metabolizer", annotationReport.getGenotypes().stream()
        .flatMap((g) -> g.getDiplotypes().stream())
        .filter((d) -> d.getGene().equals("CYP2D6"))
        .flatMap((d) -> d.getPhenotypes().stream())
        .collect(Collectors.joining("; ")));
  }

  /**
   * Can we use activity scores in outside call files? It should be specified in the column for "phenotype"
   */
  @Test
  void testOutsideCallActivityScore(TestInfo testInfo) throws Exception {
    Path outDir = TestUtils.getTestOutputDir(testInfo, false);
    Path outsideCallPath = outDir.resolve(TestUtils.getTestName(testInfo) + ".tsv");
    try (PrintWriter writer = new PrintWriter(Files.newBufferedWriter(outsideCallPath))) {
      writer.println("CYP2D6\t\t\t1.25");
    }

    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("CYP2C19");
    testWrapper.executeWithOutsideCalls(outsideCallPath);

    testWrapper.testCalledByMatcher("CYP2C19");
    testWrapper.testPrintCalls(DataSource.CPIC, "CYP2C19", "*38/*38");

    testWrapper.testNotCalledByMatcher("CYP2D6");
    testWrapper.testPrintCalls(DataSource.CPIC, "CYP2D6", "Normal Metabolizer (1.25)");

    testWrapper.testMessageCountForGene(DataSource.CPIC, "CYP2C19", 1);
    testWrapper.testGeneHasMessage(DataSource.CPIC, "CYP2C19", "reference-allele");
    testWrapper.testGeneHasMessage(DataSource.CPIC, "CYP2D6", MessageHelper.MSG_OUTSIDE_CALL);
  }

  /**
   * This test ensures that a user can specify both a diplotype AND a phenotype from an outside call. This also tests
   * to make sure the user can override the internally known phenotype with their own phenotype assignment. *2/*10 would
   * normally be a Normal Metabolizer, but this outside call overrides it as an Intermediate Metabolizer.
   */
  @Test
  void testOutsideCallActivityScoreAndPhenotype(TestInfo testInfo) throws Exception {
    Path outDir = TestUtils.getTestOutputDir(testInfo, false);
    Path outsideCallPath = outDir.resolve(TestUtils.getTestName(testInfo) + ".tsv");
    try (PrintWriter writer = new PrintWriter(Files.newBufferedWriter(outsideCallPath))) {
      writer.println("CYP2D6\t*2/*10\tIntermediate Metabolizer\t1.25");
    }

    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("CYP2C19");
    testWrapper.executeWithOutsideCalls(outsideCallPath);

    testWrapper.testCalledByMatcher("CYP2C19");
    testWrapper.testPrintCalls(DataSource.CPIC, "CYP2C19", "*38/*38");

    testWrapper.testNotCalledByMatcher("CYP2D6");
    testWrapper.testPrintCalls(DataSource.CPIC, "CYP2D6", "*2/*10");
    GeneReport cyp2d6Report = testWrapper.getContext().getGeneReport(DataSource.CPIC, "CYP2D6");
    assertNotNull(cyp2d6Report);
    assertEquals(1, cyp2d6Report.getSourceDiplotypes().size());
    assertTrue(cyp2d6Report.getSourceDiplotypes().stream().allMatch((d) -> d.getPhenotypes().contains("Intermediate Metabolizer")));

    testWrapper.testMessageCountForGene(DataSource.CPIC, "CYP2C19", 1);
    testWrapper.testGeneHasMessage(DataSource.CPIC, "CYP2C19", "reference-allele");
    testWrapper.testGeneHasMessage(DataSource.CPIC, "CYP2D6", MessageHelper.MSG_OUTSIDE_CALL);
  }

  @Test
  void testWarfarinMissingRs12777823(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("CYP2C9")
        .reference("CYP4F2")
        .reference("VKORC1")
        .missingExtraPosition("CYP2C9", "rs12777823")
    ;
    Path vcfFile = testWrapper.execute();

    testWrapper.testCalledByMatcher("CYP2C9");
    testWrapper.testReportable("CYP2C9");

    Document document = readHtmlReport(vcfFile);
    assertNotNull(document.getElementById("CYP2C9"));

    Elements warfarinCpicDips = document.select(".cpic-guideline-warfarin .rx-dip");
    assertEquals(3, warfarinCpicDips.size());
    assertEquals(warfarinCpicDips.stream()
            .map(Element::text)
            .toList(),
        List.of("CYP2C9:*1/*1", "CYP4F2:*1/*1", "VKORC1: rs9923231 reference (C)/ rs9923231 reference (C)"));

    Elements cpicWarfarinHighlightedVars = document.select(".cpic-guideline-warfarin .rx-hl-var");
    assertEquals(1, cpicWarfarinHighlightedVars.size());
    assertEquals("rs12777823:Unknown", cpicWarfarinHighlightedVars.get(0).text());
  }

  public static Document readHtmlReport(@Nullable Path file) throws IOException {
    if (file == null) {
      throw new IOException("No file specified!");
    }
    Path reporterOutput = file.getParent().resolve(BaseConfig.getBaseFilename(file) +
        BaseConfig.REPORTER_SUFFIX + ".html");
    return Jsoup.parse(reporterOutput.toFile());
  }


  /**
   * Assert single-position genes can be called outside, and double-check that a named allele gene can also be called
   * outside. Made in response to <a href="https://github.com/PharmGKB/PharmCAT/issues/154#top">issue 154</a>.
   */
  @Test
  void testOutsideSinglePositionCalls(TestInfo testInfo) throws Exception {
    Path outsideCallPath = TestUtils.createTestFile(testInfo, ".tsv");
    try (PrintWriter writer = new PrintWriter(Files.newBufferedWriter(outsideCallPath))) {
      // diplotype in backwards order
      writer.println("IFNL3\trs12979860 reference (C)/rs12979860 reference (C)\n" +
          "CYP4F2\t*1/*3");
    }

    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("CYP2C9");
    testWrapper.executeWithOutsideCalls(outsideCallPath);

    // these are outside calls
    testWrapper.testNotCalledByMatcher("IFNL3");
    testWrapper.testNotCalledByMatcher("CYP4F2");
    // this is a regular call
    testWrapper.testCalledByMatcher("CYP2C9");

    testWrapper.testPrintCpicCalls( "IFNL3", "rs12979860 reference (C)/rs12979860 reference (C)");
    testWrapper.testPrintCpicCalls( "CYP4F2", "*1/*3");
  }


  /**
   * This test ensures diplotype-specific recommendations override more generic phenotype-specific ones.
   * Specifically, phenytoin has a recommendation for "poor metabolizers" which *2/*2 will match. However, it also has a
   * recommendation for *2/*2 specifically which should "override" the poor metabolizer recommendation for this specific
   * diplotype.
   */
  @Test
  void testDiplotypeOverrideRecommendation(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("CYP2C9")
        .variation("CYP2C9", "rs1799853", "T", "T");
    testWrapper.execute();

    testWrapper.testCalledByMatcher("CYP2C9");
    testWrapper.testPrintCpicCalls( "CYP2C9", "*2/*2");
    testWrapper.testSourceDiplotypes(DataSource.CPIC, "CYP2C9", List.of("*2/*2"));
    testWrapper.testSourceDiplotypes(DataSource.DPWG, "CYP2C9", List.of("*2/*2"));

    DrugReport phenytoin = testWrapper.getContext().getDrugReport(PrescribingGuidanceSource.DPWG_GUIDELINE, "phenytoin");
    assertNotNull(phenytoin);
    List<AnnotationReport> recs = phenytoin.getGuidelines().stream().flatMap((g) -> g.getAnnotations().stream()).toList();
    assertEquals(1, recs.size());
    AnnotationReport matchingRec = recs.get(0);
    assertEquals("DPWG-PA166299254", recs.get(0).getLocalId());
  }


  @Test
  void testDuplicateEntrySecondIsBad(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .allowUnknownAllele()
        .variation("NUDT15", "rs746071566", "GAGTCG(3)", "GAGTCG(4)")
        .duplicatePositionAsIs("NUDT15", "chr13", 48037782, "0/1", "A", "C")
    ;
    Path vcfFile = testWrapper.execute();

    testWrapper.testCalledByMatcher("NUDT15");
    testWrapper.testReportable("NUDT15");

    Document document = readHtmlReport(vcfFile);
    Elements warnings = document.select("#chr13_48037782 .warningList li");
    assertNotNull(warnings);
    assertEquals(1, warnings.size());
    assertTrue(warnings.get(0).text().toLowerCase().contains("duplicate entry"));
    System.out.println(warnings.get(0).text());
  }

  @Test
  void testDuplicateEntryFirstIsBad(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .allowUnknownAllele()
        .variationAsIs("NUDT15", "chr13", 48037782, "0/1", "A", "C")
        .duplicatePosition("NUDT15", "rs746071566", "GAGTCG(3)", "GAGTCG(4)")
    ;
    Path vcfFile = testWrapper.execute();

    testWrapper.testCalledByMatcher("NUDT15");
    testWrapper.testReportable("NUDT15");

    Document document = readHtmlReport(vcfFile);
    Elements warnings = document.select("#chr13_48037782 .warningList li");
    assertNotNull(warnings);
    warnings.forEach(System.out::println);
    assertEquals(2, warnings.size());
    assertTrue(warnings.get(0).text().toLowerCase().contains("does not match expected reference"));
    assertTrue(warnings.get(1).text().toLowerCase().contains("duplicate entry"));
  }

  @Test
  void phaseSet(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false)
        .saveIntermediateFiles();

    testWrapper.getVcfBuilder()
        //chr1	97078987	.	G	T	61.6	PASS	.	GT:GQ:DP:AD:VAF:PL:PS	0|1:60:45:28,17:0.377778:61,0,66:96997594
        //chr1	97883329	.	A	G	66.2	PASS	.	GT:GQ:DP:AD:VAF:PL:PS	0|1:64:43:19,24:0.55814:66,0,67:97710720
        .reference("DPYD")
        .variationInPhaseSet("DPYD", "rs114096998", 1, "G", "T")
        .variationInPhaseSet("DPYD", "rs1801265", 1, "A", "G")
    ;

    Path vcfFile  = PathUtils.getPathToResource("org/pharmgkb/pharmcat/PipelineTest-phaseSet.vcf");
    vcfFile = testWrapper.execute(vcfFile, null, null, false);

    Document document = readHtmlReport(vcfFile);

  }
}

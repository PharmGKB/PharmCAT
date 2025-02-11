package org.pharmgkb.pharmcat;

import java.io.BufferedWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Iterator;
import java.util.List;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;
import com.google.common.collect.ImmutableSortedMap;
import com.google.common.collect.Ordering;
import org.jsoup.nodes.Document;
import org.jsoup.select.Elements;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.TestInfo;
import org.pharmgkb.pharmcat.reporter.MessageHelper;
import org.pharmgkb.pharmcat.reporter.TextConstants;
import org.pharmgkb.pharmcat.reporter.handlebars.ReportHelpers;
import org.pharmgkb.pharmcat.reporter.model.DataSource;
import org.pharmgkb.pharmcat.reporter.model.PrescribingGuidanceSource;
import org.pharmgkb.pharmcat.reporter.model.result.Diplotype;
import org.pharmgkb.pharmcat.reporter.model.result.DrugReport;
import org.pharmgkb.pharmcat.reporter.model.result.GeneReport;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.contains;
import static org.junit.jupiter.api.Assertions.*;
import static org.pharmgkb.pharmcat.PipelineTest.*;


/**
 * This JUnit test validates CYP2D6 through the pipeline.
 * This should test the data generated from a full run of the PharmCAT matcher and reporter.
 *
 * @author Mark Woon
 */
public class Cyp2d6Test {
  private static Path s_outsideCallFilePath;

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
  }

  @AfterEach
  void deleteDirectory(TestInfo testInfo) {
//    TestUtils.deleteTestOutputDirectory(testInfo);
  }


  @Test
  void testCyp2c19s4s17(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("CYP2C19", "rs12248560", "C", "T")
        .variation("CYP2C19", "rs28399504", "A", "G")
        .variation("CYP2C19", "rs3758581", "G", "G");
    testWrapper.executeWithOutsideCalls(s_outsideCallFilePath);

    testWrapper.testCalledByMatcher("CYP2C19");
    testWrapper.testReportable("CYP2D6");

    testWrapper.testPrintCpicCalls("CYP2D6", "*1/*4");
    testWrapper.testPrintCpicCalls("CYP2C19", "*1/*4");

    GeneReport cyp2d6Report = testWrapper.getContext().getGeneReport(DataSource.CPIC, "CYP2D6");
    assertNotNull(cyp2d6Report);
    assertTrue(cyp2d6Report.isOutsideCall());
  }


  @Test
  void testAmitriptylineCallWoCyp2c19(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("DPYD");
    testWrapper.executeWithOutsideCalls(s_outsideCallFilePath);

    testWrapper.testReportable("CYP2D6");
    testWrapper.testPrintCpicCalls("CYP2D6", "*1/*4");
    testWrapper.testRecommendedDiplotypes("CYP2D6", "*1", "*4");

    testWrapper.testMatchedAnnotations("amitriptyline", 3);
    testWrapper.testAnyMatchFromSource("amitriptyline", PrescribingGuidanceSource.CPIC_GUIDELINE);
    testWrapper.testAnyMatchFromSource("amitriptyline", PrescribingGuidanceSource.DPWG_GUIDELINE);
    testWrapper.testAnyMatchFromSource("amitriptyline", PrescribingGuidanceSource.FDA_ASSOC);
  }


  /**
   * Tests whether an allele that is unknown to PharmCAT/CPIC will still go through the system without throwing an
   * exception and will be reported properly.
   */
  @Test
  void testCyp2d6AlleleWithNoFunction(TestInfo testInfo) throws Exception {
    Path outsideCallPath = TestUtils.createTestFile(testInfo,".tsv");
    try (PrintWriter writer = new PrintWriter(Files.newBufferedWriter(outsideCallPath))) {
      writer.println("CYP2D6\t*1/*XXX");
    }

    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("CYP2C19");
    testWrapper.executeWithOutsideCalls(outsideCallPath);

    testWrapper.testPrintCalls(DataSource.CPIC, "CYP2D6", "*1/*XXX");
    testWrapper.testPrintCalls(DataSource.DPWG, "CYP2D6", "*1/*XXX");
    testWrapper.testSourcePhenotype(DataSource.CPIC, "CYP2D6", "n/a");

    // this nonsense allele will still match to "Indeterminate" phenotypes in guidelines for CYP2D6
    testWrapper.testMatchedAnnotations("atomoxetine", PrescribingGuidanceSource.CPIC_GUIDELINE, 2);
    DrugReport atoReport = testWrapper.getContext().getDrugReport(PrescribingGuidanceSource.CPIC_GUIDELINE, "atomoxetine");
    assertNotNull(atoReport);
    assertNotNull(atoReport.getGuidelines());
    assertEquals(1, atoReport.getGuidelines().size());
    assertTrue(atoReport.getGuidelines().stream()
        .flatMap((g) -> g.getAnnotations().stream())
        .flatMap((g) -> g.getGenotypes().stream())
        .flatMap((g) -> g.getDiplotypes().stream())
        .filter((d) -> d.getGene().equals("CYP2D6"))
        .allMatch((d) -> d.getPhenotypes().contains(TextConstants.NA)));

    GeneReport geneReport = testWrapper.getContext().getGeneReport(DataSource.CPIC, "CYP2D6");
    assertNotNull(geneReport);
    assertEquals(1, geneReport.getRecommendationDiplotypes().size());
    Diplotype diplotype = geneReport.getRecommendationDiplotypes().first();
    assertEquals("One 1.0 (Normal function) allele and one Unassigned function allele", ReportHelpers.gsFunction(diplotype));
  }


  /**
   * Should have call multimatch message.
   */
  @Test
  void testCyp2d6EquivalentDoubleCall(TestInfo testInfo) throws Exception {
    Path outsideCallPath = TestUtils.createTestFile(testInfo,".tsv");
    try (PrintWriter writer = new PrintWriter(Files.newBufferedWriter(outsideCallPath))) {
      writer.println("CYP2D6\t*1/*1");
      writer.println("CYP2D6\t*1/*2");
    }

    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("CYP2C19");
    Path vcfFile = testWrapper.executeWithOutsideCalls(outsideCallPath);

    GeneReport geneReport = testWrapper.getContext().getGeneReport(DataSource.CPIC, "CYP2D6");
    assertNotNull(geneReport);
    assertEquals(2, geneReport.getRecommendationDiplotypes().size());

    Diplotype diplotype = geneReport.getRecommendationDiplotypes().first();
    assertThat(diplotype.getPhenotypes(), contains("Normal Metabolizer"));
    assertEquals("Two 1.0 (Normal function) alleles", ReportHelpers.gsFunction(diplotype));

    Document document = readHtmlReport(vcfFile);
    Elements clomipramineSection = document.select(".guideline.clomipramine");
    assertEquals(1, clomipramineSection.size());
    assertEquals(1, clomipramineSection.get(0).getElementsByClass(MessageHelper.MSG_MULTI_CALL).size());
  }


  /**
   * Should not have call multimatch message.
   * Should have inferred CYP2D6 copy number.
   * This test has two outside calls for 2D6 with different phenotypes
   */
  @Test
  void testCyp2d6DoubleCall(TestInfo testInfo) throws Exception {
    Path outsideCallPath = TestUtils.createTestFile(testInfo,".tsv");
    try (PrintWriter writer = new PrintWriter(Files.newBufferedWriter(outsideCallPath))) {
      writer.println("CYP2D6\t*1/*1");
      writer.println("CYP2D6\t*1/*1x7");
    }

    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("CYP2C19");
    Path vcfFile = testWrapper.executeWithOutsideCalls(outsideCallPath);

    GeneReport geneReport = testWrapper.getContext().getGeneReport(DataSource.CPIC, "CYP2D6");
    assertNotNull(geneReport);
    assertEquals(2, geneReport.getRecommendationDiplotypes().size());

    Iterator<Diplotype> recommendationDipIt = geneReport.getRecommendationDiplotypes().iterator();
    Diplotype diplotype = recommendationDipIt.next();
    assertThat(diplotype.getPhenotypes(), contains("Normal Metabolizer"));
    assertEquals("Two 1.0 (Normal function) alleles", ReportHelpers.gsFunction(diplotype));

    diplotype = recommendationDipIt.next();
    assertNotNull(diplotype.getAllele2());
    assertThat(diplotype.getPhenotypes(), contains("Ultrarapid Metabolizer"));
    assertEquals("*1x" + TextConstants.GTE + "3", diplotype.getAllele2().getName());
    assertEquals("One 1.0 (Normal function) allele and one ≥3.0 (Increased function) allele", ReportHelpers.gsFunction(diplotype));

    Document document = readHtmlReport(vcfFile);
    Elements clomipramineSection = document.select(".guideline.clomipramine");
    assertEquals(1, clomipramineSection.size());
    assertEquals(0, clomipramineSection.get(0).getElementsByClass(MessageHelper.MSG_MULTI_CALL).size());
  }


  /**
   * Check that a single CYP2D6 phenotype has all the activity scores that could possibly be mapped to it.
   * Specifically, "Intermediate Metabolizer" maps to both "0.25", "0.5", "0.75" and "1.0" activity scores for CYP2D6.
   */
  @Test
  void testCyp2d6PhenotypeHasMultipleActivityScore(TestInfo testInfo) throws Exception {
    Path outDir = TestUtils.getTestOutputDir(testInfo, false);
    Path outsideCallPath = outDir.resolve(TestUtils.getTestName(testInfo) + ".tsv");
    try (PrintWriter writer = new PrintWriter(Files.newBufferedWriter(outsideCallPath))) {
      writer.println("CYP2D6\t\tIntermediate Metabolizer");
    }

    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("CYP2C19");
    Path vcfFile = testWrapper.executeWithOutsideCalls(outsideCallPath);

    List<String> cyp2c19ExpectedCalls = List.of("*38/*38");
    testWrapper.testCalledByMatcher("CYP2C19");
    testWrapper.testSourceDiplotypes(DataSource.CPIC, "CYP2C19", cyp2c19ExpectedCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "CYP2C19",
        expectedCallsToRecommendedDiplotypes(cyp2c19ExpectedCalls));
    testWrapper.testPrintCalls(DataSource.CPIC, "CYP2C19", cyp2c19ExpectedCalls);


    testWrapper.testNotCalledByMatcher("CYP2D6");
    testWrapper.testPrintCalls(DataSource.CPIC, "CYP2D6", "Intermediate Metabolizer");

    testWrapper.testMessageCountForGene(DataSource.CPIC, "CYP2C19", 1);
    testWrapper.testGeneHasMessage(DataSource.CPIC, "CYP2C19", "reference-allele");
    testWrapper.testGeneHasMessage(DataSource.CPIC, "CYP2D6", MessageHelper.MSG_OUTSIDE_CALL);

    // Check that a single CYP2D6 phenotype has all the activity scores that could possibly be mapped to it.
    // Specifically, "Intermediate Metabolizer" maps to both "0.25", "0.5", "0.75" and "1.0" activity scores for CYP2D6.
    GeneReport geneReport = testWrapper.getContext().getGeneReport(DataSource.CPIC, "CYP2D6");
    assertNotNull(geneReport);
    // we expect one diplotype to exist with multiple lookup keys
    assertEquals(1, geneReport.getRecommendationDiplotypes().size());
    Diplotype diplotype = geneReport.getRecommendationDiplotypes().iterator().next();
    // there should be no single activity score specified since this phenotype maps to more than one
    assertTrue(TextConstants.isUnspecified(diplotype.getActivityScore()));
    // there should be two and only two lookup keys, one for each activity score
    assertEquals(4, diplotype.getLookupKeys().size());
    // the two lookup keys should be the two activity scores that correspond to Intermediate Metabolizer
    assertThat(diplotype.getLookupKeys(), contains("0.25", "0.5", "0.75", "1.0"));

    Document document = readHtmlReport(vcfFile);
    SortedSet<String> expectedActivityScores = new TreeSet<>(List.of("0.25", "0.5", "0.75", "1.0"));
    htmlChecks(document,
        new ImmutableSortedMap.Builder<String, List<String>>(Ordering.natural())
            .put("CYP2C19", cyp2c19ExpectedCalls)
            .put("CYP2D6", NO_OUTSIDE_DIPLOTYPE)
            .build(),
        null, "amitriptyline", RecPresence.YES,
        new ImmutableSortedMap.Builder<String, String>(Ordering.natural())
            .put("CYP2C19", "Normal Metabolizer")
            .put("CYP2D6", "Intermediate Metabolizer")
            .build(),
        new ImmutableSortedMap.Builder<String, SortedSet<String>>(Ordering.natural())
            .put("CYP2C19", new TreeSet<>(List.of("N/A")))
            .put("CYP2D6", expectedActivityScores)
            .build(),
        RecPresence.YES,
        new ImmutableSortedMap.Builder<String, String>(Ordering.natural())
            .put("CYP2D6", "Intermediate Metabolizer")
            .build(),
        new ImmutableSortedMap.Builder<String, SortedSet<String>>(Ordering.natural())
            .put("CYP2D6", expectedActivityScores)
            .build()
    );

    assertNotNull(document.getElementById("CYP2D6"));
    Elements cyp2d6Section = document.select(".gene.CYP2D6");
    assertEquals(1, cyp2d6Section.size());
    assertEquals(1, cyp2d6Section.get(0).getElementsByClass(MessageHelper.MSG_OUTSIDE_CALL).size());
    assertEquals(0, cyp2d6Section.get(0).getElementsByClass(MessageHelper.MSG_CYP2D6_MODE).size());

    Elements clomipramineSection = document.select(".guideline.clomipramine");
    assertEquals(1, clomipramineSection.size());
    assertEquals(1, clomipramineSection.get(0).getElementsByClass(MessageHelper.MSG_MULTI_SCORE).size());
  }


  @Test
  void testCyp2d6CpicVsDpwg(TestInfo testInfo) throws Exception {
    // NOTE: this test has multiple annotations for a single population - amitriptyline
    Path outsideCallPath = TestUtils.createTestFile(testInfo,".tsv");
    try (PrintWriter writer = new PrintWriter(Files.newBufferedWriter(outsideCallPath))) {
      writer.println("CYP2D6\t*1x2/*9");
      writer.println("CYP2D6\t*1x2/*10");
      writer.println("CYP2D6\t*1x2/*17");
      writer.println("CYP2D6\t*1x3/*1");
      writer.println("CYP2D6\t*1/*1");
      writer.println("CYP2D6\t*4/*10");
      writer.println("CYP2D6\t*4/*4");
    }

    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("CYP2C19");
    testWrapper.executeWithOutsideCalls(outsideCallPath);

    List<String> expectedCyp2d6Calls = List.of("*1/*1", "*1x2/*9", "*1x2/*10", "*1x2/*17", "*1/*1x3", "*4/*4", "*4/*10");

    testWrapper.testCalledByMatcher("CYP2C19");
    testWrapper.testNotCalledByMatcher("CYP2D6");
    testWrapper.testSourceDiplotypes(DataSource.CPIC, "CYP2D6", expectedCyp2d6Calls);
    testWrapper.testPrintCalls(DataSource.CPIC, "CYP2D6", expectedCyp2d6Calls);
    testWrapper.testPrintCalls(DataSource.DPWG, "CYP2D6", expectedCyp2d6Calls);

    // TODO: finish this!
  }


  /**
   * Added to have an example of running in CYP2D6-matching mode and make sure messages are applied
   */
  @Test
  void testCyp2d6Matcher(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false, true, true);
    testWrapper.getVcfBuilder()
        .reference("CYP2D6")
        .reference("CYP2C19")
        .variation("CYP2C19", "rs3758581", "G", "G");
    Path vcfFile = testWrapper.execute();
    testWrapper.testCalledByMatcher("CYP2C19", "CYP2D6");
    testWrapper.testReportable("CYP2C19", "CYP2D6");

    Document document = readHtmlReport(vcfFile);
    assertNotNull(document.getElementById("CYP2D6"));
    Elements cyp2d6Section = document.select(".gene.CYP2D6");
    assertEquals(1, cyp2d6Section.size());
    assertEquals(0, cyp2d6Section.get(0).getElementsByClass(MessageHelper.MSG_OUTSIDE_CALL).size());
    assertEquals(1, cyp2d6Section.get(0).getElementsByClass(MessageHelper.MSG_CYP2D6_MODE).size());
  }


  @Test
  void combinationsInOutsideCalls(TestInfo testInfo) throws Exception {
    SortedMap<String, String> diplotypes = new TreeMap<>();
    diplotypes.put("*13 + *1/*36 + *10x2", "Normal Metabolizer");
    diplotypes.put("*13 + *1x2/*36 + *10", "Normal Metabolizer");
    diplotypes.put("*13 + *2/*36x2 + *10", "Normal Metabolizer");
    diplotypes.put("*13 + *68x2 + *4/*68 + *2", "Intermediate Metabolizer");
    diplotypes.put("*36x2 + *10x2/*68 + *4", "Intermediate Metabolizer");

    for (String diplotype : diplotypes.keySet()) {
      PipelineWrapper testWrapper = checkCombination(testInfo, diplotype, diplotypes.get(diplotype), true);

      if (diplotype.equals("*13 + *1/*36 + *10x2")) {
        testWrapper.testMatchedAnnotations("atomoxetine", PrescribingGuidanceSource.CPIC_GUIDELINE, 2);
        DrugReport atoReport = testWrapper.getContext().getDrugReport(PrescribingGuidanceSource.CPIC_GUIDELINE, "atomoxetine");
        assertNotNull(atoReport);
        assertNotNull(atoReport.getGuidelines());
        assertEquals(1, atoReport.getGuidelines().size());
        assertTrue(atoReport.getGuidelines().stream()
            .flatMap((g) -> g.getAnnotations().stream())
            .flatMap((g) -> g.getGenotypes().stream())
            .flatMap((g) -> g.getDiplotypes().stream())
            .filter((d) -> d.getGene().equals("CYP2D6"))
            .allMatch((d) -> d.getPhenotypes().contains("Normal Metabolizer")));

        GeneReport geneReport = testWrapper.getContext().getGeneReport(DataSource.CPIC, "CYP2D6");
        assertNotNull(geneReport);
        assertEquals(1, geneReport.getRecommendationDiplotypes().size());
        Diplotype recommendationDiplotype = geneReport.getRecommendationDiplotypes().first();
        assertEquals("One 0.5 (Decreased function) allele and one 1.0 (Normal function) allele",
            ReportHelpers.gsFunction(recommendationDiplotype));
      }
    }
  }


  PipelineWrapper checkCombination(TestInfo testInfo, String diplotype, String phenotype, boolean filenameFromDiplotype)
      throws Exception {
    Path outsideCallPath;
    String filename = null;
    if (filenameFromDiplotype) {
      filename = diplotype.replaceAll("\\*", "s")
          .replaceAll("\\[", "-")
          .replaceAll("]", "-")
          .replaceAll("/", "--")
          .replaceAll("\\+", " ")
          .replaceAll(" +", "_");
      outsideCallPath = TestUtils.createTestFile(testInfo,filename + ".tsv");
    } else {
      outsideCallPath = TestUtils.createTestFile(testInfo,".tsv");
    }
    try (PrintWriter writer = new PrintWriter(Files.newBufferedWriter(outsideCallPath))) {
      writer.println("CYP2D6\t" + diplotype);
    }

    PipelineWrapper testWrapper;
    if (filenameFromDiplotype) {
      testWrapper = new PipelineWrapper(testInfo, filename, false, false, false);
    } else {
      testWrapper = new PipelineWrapper(testInfo, false);
    }
    testWrapper.getVcfBuilder()
        .reference("CYP2C19");
    testWrapper.executeWithOutsideCalls(outsideCallPath);

    testWrapper.testPrintCalls(DataSource.CPIC, "CYP2D6", diplotype);
    testWrapper.testSourcePhenotype(DataSource.CPIC, "CYP2D6", phenotype);
    return testWrapper;
  }
}

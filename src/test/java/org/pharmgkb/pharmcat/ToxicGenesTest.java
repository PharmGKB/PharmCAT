package org.pharmgkb.pharmcat;

import java.nio.file.Path;
import java.util.List;
import java.util.SortedMap;
import java.util.TreeMap;
import com.google.common.collect.ImmutableSortedMap;
import com.google.common.collect.Ordering;
import org.jsoup.nodes.Document;
import org.jsoup.select.Elements;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.TestInfo;
import org.pharmgkb.pharmcat.reporter.TextConstants;
import org.pharmgkb.pharmcat.reporter.handlebars.ReportHelpers;
import org.pharmgkb.pharmcat.reporter.model.DataSource;
import org.pharmgkb.pharmcat.reporter.model.PrescribingGuidanceSource;
import org.pharmgkb.pharmcat.reporter.model.result.DrugReport;
import org.pharmgkb.pharmcat.reporter.model.result.GeneReport;
import org.pharmgkb.pharmcat.reporter.model.result.Genotype;
import org.pharmgkb.pharmcat.reporter.model.result.GuidelineReport;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.pharmgkb.pharmcat.PipelineTest.*;


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
  void testUncallable_partial_haplotype(TestInfo testInfo) throws Exception {

    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("TPMT", "rs1256618794", "C", "A") // C -> A
        .variation("TPMT", "rs753545734", "C", "C") // C -> T
    ;
    Path vcfFile = testWrapper.execute();

    List<String> expectedCalls = UNKNOWN_CALL;

    testWrapper.testNotCalledByMatcher("TPMT");
    testWrapper.testSourceDiplotypes(DataSource.CPIC, "TPMT", expectedCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "TPMT", expectedCalls);
    testWrapper.testPrintCalls(DataSource.CPIC, "TPMT", List.of(TextConstants.UNCALLED));

    Document document = readHtmlReport(vcfFile);
    SortedMap<String, List<String>> expectedCallsMap = new TreeMap<>();
    expectedCallsMap.put("TPMT", UNKNOWN_CALL);
    htmlCheckGenes(document,
        new ImmutableSortedMap.Builder<String, List<String>>(Ordering.natural())
            .put("TPMT", UNKNOWN_CALL)
            .build(),
        null);
  }


  @Test
  void cacna1sRef_ryr1Ref(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("CACNA1S")
        .reference("RYR1");
    Path vcfFile = testWrapper.execute();

    List<String> expectedCalls = List.of(TextConstants.HOMOZYGOUS_REFERENCE);

    testWrapper.testCalledByMatcher("CACNA1S", "RYR1");

    testWrapper.testSourceDiplotypes(DataSource.CPIC, "CACNA1S", expectedCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "CACNA1S", expectedCallsToRecommendedDiplotypes(expectedCalls));
    testWrapper.testPrintCalls(DataSource.CPIC, "CACNA1S", expectedCalls);

    testWrapper.testSourceDiplotypes(DataSource.CPIC, "RYR1", expectedCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, "RYR1", expectedCallsToRecommendedDiplotypes(expectedCalls));
    testWrapper.testPrintCalls(DataSource.CPIC, "RYR1", expectedCalls);

    // each gene has its own annotation which means there should be 2 CPIC annotations matched, one for each gene
    testWrapper.testMatchedAnnotations("desflurane", PrescribingGuidanceSource.CPIC_GUIDELINE, 1);
    testWrapper.testNoMatchFromSource("desflurane", PrescribingGuidanceSource.DPWG_GUIDELINE);

    Document document = readHtmlReport(vcfFile);
    htmlChecks(document, new ImmutableSortedMap.Builder<String, List<String>>(Ordering.natural())
            .put("CACNA1S", expectedCalls)
            .put("RYR1", expectedCalls)
            .build(),
        "desflurane", RecPresence.YES, RecPresence.NO);
  }


  @Test
  void testG6pdRef_male(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .male()
        .reference("G6PD");
    Path vcfFile = testWrapper.execute();

    testWrapper.testCalledByMatcher("G6PD");
    testWrapper.testReportable("G6PD");

    Document document = readHtmlReport(vcfFile);
    Elements g6pdSections = document.select(".gene.G6PD");
    assertEquals(1, g6pdSections.size());
    Elements g6pdCallElems = g6pdSections.get(0).getElementsByClass("genotype-result");
    assertEquals(1, g6pdCallElems.size());
    assertEquals("B (reference)", g6pdCallElems.text());
  }

  @Test
  void testG6pd_Ref_female(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .female()
        .reference("G6PD");
    Path vcfFile = testWrapper.execute();

    testWrapper.testCalledByMatcher("G6PD");
    testWrapper.testReportable("G6PD");

    Document document = readHtmlReport(vcfFile);
    Elements g6pdSections = document.select(".gene.G6PD");
    assertEquals(1, g6pdSections.size());
    Elements g6pdCallElems = g6pdSections.get(0).getElementsByClass("genotype-result");
    assertEquals(1, g6pdCallElems.size());
    assertEquals("B (reference)/B (reference)", g6pdCallElems.text());
  }

  @Test
  void testG6pd_MDPSCB_male(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .male()
        .reference("G6PD")
        .variation("G6PD", "rs5030868", "A");
    Path vcfFile = testWrapper.execute();

    String gene = "G6PD";
    List<String> expectedCalls = List.of("Mediterranean, Dallas, Panama, Sassari, Cagliari, Birmingham");

    testWrapper.testCalledByMatcher(gene);
    testWrapper.testSourceDiplotypes(DataSource.CPIC, gene, expectedCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, gene, expectedCalls);
    testWrapper.testReportable(gene);
    testWrapper.testPrintCalls(DataSource.CPIC, gene, expectedCalls);

    Document document = readHtmlReport(vcfFile);
    htmlChecks(document, gene, expectedCalls, "aspirin", RecPresence.NO, RecPresence.NO);
  }

  @Test
  void testG6pd_MDPSCB_female_homo(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("G6PD")
        .variation("G6PD", "rs5030868", "A", "A");
    Path vcfFile = testWrapper.execute();

    String gene = "G6PD";
    List<String> expectedCalls = List.of("Mediterranean, Dallas, Panama, Sassari, Cagliari, Birmingham/Mediterranean, Dallas, Panama, Sassari, Cagliari, Birmingham");

    testWrapper.testCalledByMatcher(gene);
    testWrapper.testSourceDiplotypes(DataSource.CPIC, gene, expectedCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, gene, expectedCallsToRecommendedDiplotypes(expectedCalls));
    testWrapper.testReportable(gene);
    testWrapper.testPrintCalls(DataSource.CPIC, gene, expectedCalls);

    Document document = readHtmlReport(vcfFile);
    htmlChecks(document, gene, expectedCalls, "aspirin", RecPresence.NO, RecPresence.NO);
  }

  @Test
  void testG6pd_MDPSCB_female_het(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("G6PD")
        .variation("G6PD", "rs5030868", "G", "A");
    Path vcfFile = testWrapper.execute();

    String gene = "G6PD";
    List<String> expectedCalls = List.of("B (reference)/Mediterranean, Dallas, Panama, Sassari, Cagliari, Birmingham");

    testWrapper.testCalledByMatcher(gene);
    testWrapper.testSourceDiplotypes(DataSource.CPIC, gene, expectedCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, gene, expectedCallsToRecommendedDiplotypes(expectedCalls));
    testWrapper.testReportable(gene);
    testWrapper.testPrintCalls(DataSource.CPIC, gene, expectedCalls);

    Document document = readHtmlReport(vcfFile);
    htmlChecks(document, gene, expectedCalls, "aspirin", RecPresence.NO, RecPresence.NO);
  }


  @Test
  void testG6pd_chatham_male(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .male()
        .reference("G6PD")
        .variation("G6PD", "rs5030869", "T");
    Path vcfFile = testWrapper.execute();

    String gene = "G6PD";
    List<String> expectedCalls = List.of("Chatham");

    testWrapper.testCalledByMatcher(gene);
    testWrapper.testSourceDiplotypes(DataSource.CPIC, gene, expectedCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, gene, expectedCalls);
    testWrapper.testReportable(gene);
    testWrapper.testPrintCalls(DataSource.CPIC, gene, expectedCalls);

    Document document = readHtmlReport(vcfFile);
    htmlChecks(document, gene, expectedCalls, "aspirin", RecPresence.NO, RecPresence.NO);
  }

  @Test
  void testG6pd_chatham_female_homo(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("G6PD")
        .variation("G6PD", "rs5030869", "T", "T");
    Path vcfFile = testWrapper.execute();

    String gene = "G6PD";
    List<String> expectedCalls = List.of("Chatham/Chatham");

    testWrapper.testCalledByMatcher(gene);
    testWrapper.testSourceDiplotypes(DataSource.CPIC, gene, expectedCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, gene, expectedCallsToRecommendedDiplotypes(expectedCalls));
    testWrapper.testReportable(gene);
    testWrapper.testPrintCalls(DataSource.CPIC, gene, expectedCalls);

    Document document = readHtmlReport(vcfFile);
    htmlChecks(document, gene, expectedCalls, "aspirin", RecPresence.NO, RecPresence.NO);
  }

  @Test
  void testG6pd_chatham_female_het(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("G6PD")
        .variation("G6PD", "rs5030869", "C", "T");
    Path vcfFile = testWrapper.execute();

    String gene = "G6PD";
    List<String> expectedCalls = List.of("B (reference)/Chatham");

    testWrapper.testCalledByMatcher(gene);
    testWrapper.testSourceDiplotypes(DataSource.CPIC, gene, expectedCalls);
    testWrapper.testRecommendedDiplotypes(DataSource.CPIC, gene, expectedCallsToRecommendedDiplotypes(expectedCalls));
    testWrapper.testReportable(gene);
    testWrapper.testPrintCalls(DataSource.CPIC, gene, expectedCalls);

    Document document = readHtmlReport(vcfFile);
    htmlChecks(document, gene, expectedCalls, "aspirin", RecPresence.NO, RecPresence.NO);
  }

  @Test
  void testG6pd_Arakawa_male(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .male()
        .reference("G6PD")
        .variation("G6PD", "chrX", 154532082, "A");
    Path vcfFile = testWrapper.execute();

    testWrapper.testCalledByMatcher("G6PD");
    testWrapper.testReportable("G6PD");

    Document document = readHtmlReport(vcfFile);
    Elements g6pdSections = document.select(".gene.G6PD");
    assertEquals(1, g6pdSections.size());
    Elements g6pdCallElems = g6pdSections.get(0).getElementsByClass("genotype-result");
    assertEquals(1, g6pdCallElems.size());
    assertEquals("Arakawa", g6pdCallElems.text());
  }

  @Test
  void testG6pd_Arakawa_female_het(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .female()
        .reference("G6PD")
        .variation("G6PD", "chrX", 154532082, "G", "A");
    Path vcfFile = testWrapper.execute();

    testWrapper.testCalledByMatcher("G6PD");
    testWrapper.testReportable("G6PD");

    Document document = readHtmlReport(vcfFile);
    Elements g6pdSections = document.select(".gene.G6PD");
    assertEquals(1, g6pdSections.size());
    Elements g6pdCallElems = g6pdSections.get(0).getElementsByClass("genotype-result");
    assertEquals(1, g6pdCallElems.size());
    assertEquals("Arakawa/B (reference)", g6pdCallElems.text());
  }

  @Test
  void testG6pd_Arakawa_female_homo(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .female()
        .reference("G6PD")
        .variation("G6PD", "chrX", 154532082, "A", "A");
    Path vcfFile = testWrapper.execute();

    testWrapper.testCalledByMatcher("G6PD");
    testWrapper.testReportable("G6PD");

    Document document = readHtmlReport(vcfFile);
    Elements g6pdSections = document.select(".gene.G6PD");
    assertEquals(1, g6pdSections.size());
    Elements g6pdCallElems = g6pdSections.get(0).getElementsByClass("genotype-result");
    assertEquals(1, g6pdCallElems.size());
    assertEquals("Arakawa/Arakawa", g6pdCallElems.text());
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
}

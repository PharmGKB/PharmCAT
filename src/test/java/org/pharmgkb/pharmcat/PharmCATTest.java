package org.pharmgkb.pharmcat;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Optional;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.stream.Collectors;
import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import org.jsoup.Jsoup;
import org.jsoup.nodes.Document;
import org.jsoup.nodes.Element;
import org.jsoup.select.Elements;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.TestInfo;
import org.pharmgkb.common.util.PathUtils;
import org.pharmgkb.pharmcat.haplotype.ResultSerializer;
import org.pharmgkb.pharmcat.haplotype.model.GeneCall;
import org.pharmgkb.pharmcat.haplotype.model.Result;
import org.pharmgkb.pharmcat.phenotype.Phenotyper;
import org.pharmgkb.pharmcat.reporter.BadOutsideCallException;
import org.pharmgkb.pharmcat.reporter.MessageHelper;
import org.pharmgkb.pharmcat.reporter.ReportContext;
import org.pharmgkb.pharmcat.reporter.model.DataSource;
import org.pharmgkb.pharmcat.reporter.model.MessageAnnotation;
import org.pharmgkb.pharmcat.reporter.model.VariantReport;
import org.pharmgkb.pharmcat.reporter.model.result.AnnotationReport;
import org.pharmgkb.pharmcat.reporter.model.result.Diplotype;
import org.pharmgkb.pharmcat.reporter.model.result.DrugReport;
import org.pharmgkb.pharmcat.reporter.model.result.GeneReport;
import org.pharmgkb.pharmcat.reporter.model.result.Genotype;
import org.pharmgkb.pharmcat.reporter.model.result.GuidelineReport;
import org.pharmgkb.pharmcat.reporter.model.result.Haplotype;

import static com.github.stefanbirkner.systemlambda.SystemLambda.tapSystemOut;
import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.contains;
import static org.junit.jupiter.api.Assertions.*;


/**
 * Test the data generated from a full run of the PharmCAT matcher and reporter.
 *
 * @author Ryan Whaley
 */
class PharmCATTest {
  private static final String sf_outsideCalls = """
      ##Test Outside Call Data
      #Gene\tDiplotype\tdiplotype activity\tdiplotype calling notes\tjaccard\tpart\tpValue\tROI notes\tspecial case\tnomenclature version
      CYP2D6\tCYP2D6*1/CYP2D6*4\t\t\t0.6\t0.75\tp: 0.0\t\t\tv1.9-2017_02_09
      """;
  private static final String sf_otherOutsideCalls = "CYP2D6\t*3/*4\nG6PD\tB (wildtype)/B (wildtype)\n";
  private static final String sf_diplotypesTemplate = "\nmatcher: %s\nreporter: %s\nprint (displayCalls): %s";
  private static boolean m_compact;
  private static List<DataSource> m_sources = Lists.newArrayList(DataSource.CPIC, DataSource.DPWG);
  private static Path s_outsideCallFilePath;
  private static Path s_otherOutsideCallFilePath;


  static void setCompact(boolean compact) {
    m_compact = compact;
  }

  static void setSources(List<DataSource> sources) {
    m_sources = sources;
  }


  @BeforeAll
  static void prepare() throws IOException {

    s_outsideCallFilePath = TestUtils.createTestFile(PharmCATTest.class, "outsideCall.tsv");
    try (FileWriter fw = new FileWriter(s_outsideCallFilePath.toFile())) {
      fw.write(sf_outsideCalls);
    }

    s_otherOutsideCallFilePath = TestUtils.createTestFile(PharmCATTest.class, "otherOutsideCall.tsv");
    try (FileWriter fw = new FileWriter(s_otherOutsideCallFilePath.toFile())) {
      fw.write(sf_otherOutsideCalls);
    }
  }


  @Test
  void cliBasic() throws Exception {
    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/reference.vcf");
    Path refMatcherOutput = vcfFile.getParent().resolve("reference.match.json");
    Path refPhenotyperOutput = vcfFile.getParent().resolve("reference.phenotype.json");
    Path refReporterOutput = vcfFile.getParent().resolve("reference.report.html");

    // require VCF
    String systemOut = tapSystemOut(() -> PharmCAT.main(null));
    //System.out.println(systemOut);
    assertTrue(systemOut.contains("No input"));
    assertTrue(systemOut.contains("-vcf"));


    // standard full run, should output to same directory as VCF file
    try {
      systemOut = tapSystemOut(() -> PharmCAT.main(new String[] {
          "-vcf", vcfFile.toString(),
      }));
      System.out.println(systemOut);
      assertTrue(systemOut.contains("Done."));
      assertTrue(Files.exists(refMatcherOutput));
      assertTrue(Files.exists(refPhenotyperOutput));
      assertTrue(Files.exists(refReporterOutput));

      ResultSerializer resultSerializer = new ResultSerializer();
      Result result = resultSerializer.fromJson(refMatcherOutput);
      Optional<GeneCall> gcOpt = result.getGeneCalls().stream()
          .filter(gc -> gc.getGene().equals("CYP2D6"))
          .findFirst();
      assertTrue(gcOpt.isEmpty());

      Collection<GeneReport> reports = Phenotyper.read(refPhenotyperOutput).getGeneReports().get(DataSource.CPIC)
          .values();
      Optional<GeneReport> grOpt = reports.stream()
          .filter(gr -> gr.getGene().equals("CYP2D6"))
          .findFirst();
      assertTrue(grOpt.isPresent());
      assertFalse(grOpt.get().isCalled());

      Document document = Jsoup.parse(refReporterOutput.toFile());
      assertNotNull(document.getElementById("section-i"));
      assertNotNull(document.getElementById("CYP2D6"));
      assertNotNull(document.getElementById("aripiprazole"));

    } finally {
      // always delete
      boolean save = TestUtils.isSaveTestOutput();
      TestUtils.setSaveTestOutput(false);
      TestUtils.deleteTestFiles(refMatcherOutput, refPhenotyperOutput, refReporterOutput);
      TestUtils.setSaveTestOutput(save);
    }
  }


  @Test
  void cliReference(TestInfo testInfo) throws Exception {
    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/reference.vcf");

    Path outputDir = TestUtils.getTestOutputDir(testInfo, true);
    Path refMatcherOutput = outputDir.resolve("reference.match.json");
    Path refPhenotyperOutput = outputDir.resolve("reference.phenotype.json");
    Path refReporterOutput = outputDir.resolve("reference.report.html");

    try {
      String systemOut = tapSystemOut(() -> PharmCAT.main(new String[] {
          "-vcf", vcfFile.toString(),
          "-o", outputDir.toString(),
      }));
      //System.out.println(systemOut);
      assertTrue(systemOut.contains("Done."));
      assertTrue(Files.exists(refMatcherOutput));
      assertTrue(Files.exists(refPhenotyperOutput));
      assertTrue(Files.exists(refReporterOutput));

      ResultSerializer resultSerializer = new ResultSerializer();
      Result result = resultSerializer.fromJson(refMatcherOutput);
      Optional<GeneCall> gcOpt = result.getGeneCalls().stream()
          .filter(gc -> gc.getGene().equals("CYP2D6"))
          .findFirst();
      assertTrue(gcOpt.isEmpty());

      Collection<GeneReport> reports = Phenotyper.read(refPhenotyperOutput).getGeneReports().get(DataSource.CPIC)
          .values();
      Optional<GeneReport> grOpt = reports.stream()
          .filter(gr -> gr.getGene().equals("CYP2D6"))
          .findFirst();
      assertTrue(grOpt.isPresent());
      assertFalse(grOpt.get().isCalled());

      Document document = Jsoup.parse(refReporterOutput.toFile());
      assertNotNull(document.getElementById("section-i"));
      assertNotNull(document.getElementById("CYP2D6"));
      assertNotNull(document.getElementById("aripiprazole"));

    } finally {
      TestUtils.deleteTestFiles(refMatcherOutput, refPhenotyperOutput, refReporterOutput);
    }
  }


  @Test
  void cliCallCyp2d6(TestInfo testInfo) throws Exception {
    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/reference.vcf");

    Path outputDir = TestUtils.getTestOutputDir(testInfo, true);
    Path refMatcherOutput = outputDir.resolve("reference.match.json");
    Path refPhenotyperOutput = outputDir.resolve("reference.phenotype.json");
    Path refReporterOutput = outputDir.resolve("reference.report.html");

    try {
      String systemOut = tapSystemOut(() -> PharmCAT.main(new String[]{
          "-vcf", vcfFile.toString(),
          "-o", outputDir.toString(),
          "-research", "cyp2d6"
      }));
      //System.out.println(systemOut);
      assertTrue(systemOut.contains("Done."));
      assertTrue(Files.exists(refMatcherOutput));
      assertTrue(Files.exists(refPhenotyperOutput));
      assertTrue(Files.exists(refReporterOutput));

      ResultSerializer resultSerializer = new ResultSerializer();
      Result result = resultSerializer.fromJson(refMatcherOutput);
      Optional<GeneCall> opt = result.getGeneCalls().stream()
          .filter(gc -> gc.getGene().equals("CYP2D6"))
          .findFirst();
      assertTrue(opt.isPresent());
      assertEquals(1, opt.get().getDiplotypes().size());
      assertEquals("*1/*1", opt.get().getDiplotypes().iterator().next().getName());

      Collection<GeneReport> reports = Phenotyper.read(refPhenotyperOutput).getGeneReports().get(DataSource.CPIC)
          .values();
      Optional<GeneReport> grOpt = reports.stream()
          .filter(gr -> gr.getGene().equals("CYP2D6"))
          .findFirst();
      assertTrue(grOpt.isPresent());
      assertTrue(grOpt.get().isCalled());
      assertFalse(grOpt.get().isOutsideCall());

      Document document = Jsoup.parse(refReporterOutput.toFile());
      Element geneTitle = document.getElementById("CYP2D6");
      assertNotNull(geneTitle);
      assertNotNull(geneTitle.parent());
      Elements diplotypes = geneTitle.parent().getElementsByTag("li");
      for (String diplotype : diplotypes.eachText()) {
        assertEquals("*1/*1", diplotype);
      }
      assertNotNull(document.getElementById("aripiprazole"));

    } finally {
      TestUtils.deleteTestFiles(outputDir);
    }
  }

  @Test
  void cliCallCyp2d6WithOverlappingOutsideCall(TestInfo testInfo) throws Exception {
    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/reference.vcf");
    Path outsideCallFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/PharmCATTest-cyp2d6.tsv");

    Path outputDir = TestUtils.getTestOutputDir(testInfo, true);
    Path refMatcherOutput = outputDir.resolve("reference.match.json");
    Path refPhenotyperOutput = outputDir.resolve("reference.phenotype.json");
    Path refReporterOutput = outputDir.resolve("reference.report.html");

    try {
      String systemOut = tapSystemOut(() -> PharmCAT.main(new String[]{
          "-vcf", vcfFile.toString(),
          "-po", outsideCallFile.toString(),
          "-o", outputDir.toString(),
          "-matcher", "-research", "cyp2d6"
      }));
      System.out.println(systemOut);
      assertTrue(Files.exists(refMatcherOutput));
      assertFalse(Files.exists(refPhenotyperOutput));
      assertFalse(Files.exists(refReporterOutput));

    } finally {
      TestUtils.deleteTestFiles(outputDir);
    }
  }

  @Test
  void cliOutsideCallsOnly(TestInfo testInfo) throws Exception {
    Path outsideCallFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/PharmCATTest-cyp2d6.tsv");

    Path outputDir = TestUtils.getTestOutputDir(testInfo, true);
    Path matcherOutput = outputDir.resolve("PharmCATTest-cyp2d6.match.json");
    Path phenotyperOutput = outputDir.resolve("PharmCATTest-cyp2d6.phenotype.json");
    Path reporterOutput = outputDir.resolve("PharmCATTest-cyp2d6.report.html");

    try {
      String systemOut = tapSystemOut(() -> PharmCAT.main(new String[] {
          "-o", outputDir.toString(),
          "-phenotyper",
          "-po", outsideCallFile.toString()
      }));
      System.out.println(systemOut);
      assertTrue(systemOut.contains("Done."));
      assertFalse(Files.exists(matcherOutput));
      assertTrue(Files.exists(phenotyperOutput));
      assertFalse(Files.exists(reporterOutput));

      Collection<GeneReport> reports = Phenotyper.read(phenotyperOutput).getGeneReports().get(DataSource.CPIC)
          .values();
      Optional<GeneReport> grOpt = reports.stream()
          .filter(gr -> gr.getGene().equals("CYP2D6"))
          .findFirst();
      assertTrue(grOpt.isPresent());
      assertFalse(grOpt.get().isCalled());
      assertTrue(grOpt.get().isOutsideCall());

    } finally {
      TestUtils.deleteTestFiles(outputDir);
    }
  }

  /**
   * This test is dependent on CYP2C19 definition.
   * It may fail if CYP2C19 definition changes too much.
   */
  @Test
  void cliMatchAllFlag(TestInfo testInfo) throws Exception {
    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/PharmCATTest-cyp2c19MissingPositions.vcf");

    Path outputDir = TestUtils.getTestOutputDir(testInfo, true);
    String baseFilename = TestUtils.getTestName(testInfo);
    Path matcherOutput = outputDir.resolve(baseFilename + ".match.json");
    Path phenotyperOutput = outputDir.resolve(baseFilename + ".phenotype.json");
    Path reporterOutput = outputDir.resolve(baseFilename + ".report.html");

    // matcher only, expecting 1 CYP2C19 matches
    try {
      String systemOut = tapSystemOut(() -> PharmCAT.main(new String[] {
          "-vcf", vcfFile.toString(),
          "-matcher",
          "-o", outputDir.toString(),
          "-bf", baseFilename,
      }));
      //System.out.println(systemOut);
      assertTrue(systemOut.contains("Done."));
      assertTrue(Files.exists(matcherOutput));
      assertFalse(Files.exists(phenotyperOutput));
      assertFalse(Files.exists(reporterOutput));

      ResultSerializer resultSerializer = new ResultSerializer();
      Result result = resultSerializer.fromJson(matcherOutput);
      Optional<GeneCall> gcOpt = result.getGeneCalls().stream()
          .filter(gc -> gc.getGene().equals("CYP2C19"))
          .findFirst();
      assertTrue(gcOpt.isPresent());
      GeneCall gc = gcOpt.get();
      assertEquals(1, gc.getDiplotypes().size());

    } finally {
      TestUtils.deleteTestFiles(outputDir);
    }

    // matcher only, expecting many CYP2C19 matches
    try {
      String systemOut = tapSystemOut(() -> PharmCAT.main(new String[] {
          "-vcf", vcfFile.toString(),
          "-matcher",
          "-ma",
          "-o", outputDir.toString(),
          "-bf", baseFilename,
      }));
      //System.out.println(systemOut);
      assertTrue(systemOut.contains("Done."));
      assertTrue(Files.exists(matcherOutput));
      assertFalse(Files.exists(phenotyperOutput));
      assertFalse(Files.exists(reporterOutput));

      ResultSerializer resultSerializer = new ResultSerializer();
      Result result = resultSerializer.fromJson(matcherOutput);
      Optional<GeneCall> gcOpt = result.getGeneCalls().stream()
          .filter(gc -> gc.getGene().equals("CYP2C19"))
          .findFirst();
      assertTrue(gcOpt.isPresent());
      GeneCall gc = gcOpt.get();
      assertTrue(gc.getDiplotypes().size() > 50);

    } finally {
      TestUtils.deleteTestFiles(outputDir);
    }
  }


  /**
   * This test ensures that the output of the phenotyper is consistent regardless of whether you run the matcher first
   * in one process then the phenotyper in a second process versus if you run the matcher and phenotyper together in
   * one process.
   */
  @Test
  void cliTestConsistentOutput(TestInfo testInfo) throws Exception {
    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/reference.vcf");

    Path outputDir = TestUtils.getTestOutputDir(testInfo, true);
    Path singlesMatcherOutput = outputDir.resolve("singles.match.json");
    Path singlesPhenotyperOutput = outputDir.resolve("singles.phenotype.json");
    Path doublePhenotyperOutput = outputDir.resolve("double.phenotype.json");

    try {
      String systemOut = tapSystemOut(() -> PharmCAT.main(new String[] {
          "-matcher",
          "-vcf", vcfFile.toString(),
          "-o", outputDir.toString(),
          "-bf", "singles",
      }));
      assertTrue(systemOut.contains("Done."));
      assertTrue(Files.exists(singlesMatcherOutput));

      String singlesPhenoOut = tapSystemOut(() -> PharmCAT.main(new String[]{
          "-phenotyper",
          "-pi", singlesMatcherOutput.toString(),
          "-o", outputDir.toString(),
          "-bf", "singles"
      }));
      assertTrue(singlesPhenoOut.contains("Done."));
      assertTrue(Files.exists(singlesPhenotyperOutput));

      String doubleOut = tapSystemOut(() -> PharmCAT.main(new String[] {
          "-matcher",
          "-phenotyper",
          "-vcf", vcfFile.toString(),
          "-o", outputDir.toString(),
          "-bf", "double",
      }));
      assertTrue(doubleOut.contains("Done."));
      assertTrue(Files.exists(doublePhenotyperOutput));

      assertEquals(Files.readString(singlesPhenotyperOutput), Files.readString(doublePhenotyperOutput));
    } finally {
      TestUtils.deleteTestFiles(outputDir);
    }
  }


  /**
   * NOTE: if these assertions fail then new data may have been added from the DataManager because of an update to the
   * CPIC database. If that's true, then update these numbers to the current count. If the count changes with no known
   * change to the CPIC database then something may be wrong in code.
   */
  @Test
  void testCounts(TestInfo testInfo) throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
            .reference("CYP2C9");
    testWrapper.execute(null);
    assertEquals(23, testWrapper.getContext().getGeneReports().keySet().stream()
        .flatMap((k) -> testWrapper.getContext().getGeneReports().get(k).values().stream()
            .map(GeneReport::getGeneDisplay))
        .collect(Collectors.toSet())
        .size()
    );
    assertEquals(125, testWrapper.getContext().getDrugReports().keySet().stream()
        .flatMap((k) -> testWrapper.getContext().getDrugReports().get(k).values().stream()
            .map(DrugReport::getName))
        .collect(Collectors.toSet())
        .size());
  }

  @Test
  void testAll(TestInfo testInfo) throws Exception {
    Path outsideCallPath = TestUtils.createTestFile(testInfo, ".tsv");
    try (PrintWriter writer = new PrintWriter(Files.newBufferedWriter(outsideCallPath))) {
      writer.println(
          """
              CYP2D6\t*3/*4
              G6PD\tB (wildtype)/B (wildtype)
              HLA-A\t\t*31:01 positive
              HLA-B\t*15:02/*57:01"""
      );
    }
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper(testInfo, false);
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
        .reference("IFNL3")
        .reference("NUDT15")
        .reference("RYR1")
        .reference("SLCO1B1")
        .reference("TPMT")
        .reference("UGT1A1")
        .reference("VKORC1");
    testWrapper.execute(outsideCallPath);

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
        "IFNL3",
        "NUDT15",
        "RYR1",
        "SLCO1B1",
        "TPMT",
        "UGT1A1",
        "VKORC1"
    );
    testWrapper.testNotCalledByMatcher("CYP2D6", "G6PD", "HLA-A", "HLA-B");
  }

  /**
   * This test illustrates when one gene in a two-gene guideline (amitriptyline) is not called that it should still be
   * able to come up with a matched annotation.
   */
  @Test
  void testCyp2c19(TestInfo testInfo) throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("CYP2C19", "rs3758581", "G", "G");
    testWrapper.execute(s_otherOutsideCallFilePath);

    testWrapper.testCalledByMatcher("CYP2C19");
    testWrapper.testPrintCpicCalls( "CYP2C19", "*1/*1");

    testWrapper.testMatchedAnnotations("amitriptyline", DataSource.CPIC, 1);
    testWrapper.testMatchedAnnotations("amitriptyline", DataSource.DPWG, 1);
    testWrapper.testMatchedAnnotations("citalopram", DataSource.CPIC, 1);
    testWrapper.testMatchedAnnotations("citalopram", DataSource.DPWG, 1);
    testWrapper.testMatchedAnnotations("ivacaftor", 0);
  }

  /**
   * Tests how PharmCAT handles that state when sample VCF data exists for a gene and an outside call also exists for
   * that gene. Currently, this should execute successfully by ignoring VCF data and using the outside call
   */
  @Test
  void testCallerCollision(TestInfo testInfo) throws Exception {
    Path outsideCallPath = TestUtils.createTestFile(testInfo, ".tsv");
    try (PrintWriter writer = new PrintWriter(Files.newBufferedWriter(outsideCallPath))) {
      writer.println("CYP2C19\t*2/*2");
    }

    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("CYP2C19");
    testWrapper.execute(outsideCallPath);

    testWrapper.testNotCalledByMatcher("CYP2C19");
    // this is the diplotype indicated in the outside call, not the one matched
    testWrapper.testPrintCpicCalls( "CYP2C19", "*2/*2");

    testWrapper.testMessageCountForGene(DataSource.CPIC, "CYP2C19", 2);
    testWrapper.testGeneHasMessage(DataSource.CPIC, "CYP2C19", "prefer-sample-data",
        MessageHelper.MSG_OUTSIDE_CALL);
  }


  /**
   * This test case demos that an "ambiguity" {@link MessageAnnotation} which specifies a variant and a diplotype call
   * for a given drug report will be matched and added to the {@link DrugReport}
   */
  @Test
  void testCyp2c19_s1s2rs58973490het(TestInfo testInfo) throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("CYP2C19", "rs12769205", "A", "G")
        .variation("CYP2C19", "rs58973490", "G", "A")
        .variation("CYP2C19", "rs4244285", "G", "A")
        .variation("CYP2C19", "rs3758581", "G", "G");
    testWrapper.execute(s_otherOutsideCallFilePath);

    testWrapper.testCalledByMatcher("CYP2C19");
    testWrapper.testPrintCpicCalls( "CYP2C19", "*1/*2");
    testWrapper.testNotCalledByMatcher("CYP2D6");
    testWrapper.testPrintCpicCalls( "CYP2D6", "*3/*4");

    testWrapper.testMatchedAnnotations("amitriptyline", 2);
    testWrapper.testAnyMatchFromSource("amitriptyline", DataSource.CPIC);
    testWrapper.testAnyMatchFromSource("amitriptyline", DataSource.DPWG);
    testWrapper.testMatchedAnnotations("citalopram", 2);
    testWrapper.testAnyMatchFromSource("citalopram", DataSource.CPIC);
    testWrapper.testAnyMatchFromSource("citalopram", DataSource.DPWG);
    testWrapper.testMatchedAnnotations("clomipramine", 3);
    testWrapper.testAnyMatchFromSource("clomipramine", DataSource.CPIC);
    testWrapper.testAnyMatchFromSource("clomipramine", DataSource.DPWG);
    testWrapper.testMatchedAnnotations("ivacaftor", 0);

    GeneReport cyp2c19report = testWrapper.getContext().getGeneReport(DataSource.CPIC, "CYP2C19");

    VariantReport vr = cyp2c19report.findVariantReport("rs58973490")
        .orElseThrow(() -> new RuntimeException("Variant missing from test data"));
    assertTrue(vr.isHetCall());

    // ambiguity message will not apply in this case because all variants are available for CYP2C19, but one message
    // should appear for the *1 call
    assertEquals(1, cyp2c19report.getMessages().stream()
        .filter(m -> m.getExceptionType().equals(MessageAnnotation.TYPE_AMBIGUITY) &&
            Objects.requireNonNull(m.getMatches().getVariant()).equals("rs58973490"))
        .count());

    testWrapper.testMessageCountForDrug(DataSource.CPIC, "amitriptyline", 1);
  }

  /**
   * This test case demos that an "ambiguity" {@link MessageAnnotation} which specifies a variant and a diplotype call
   * for a given drug report will not be matched when the variant in the message is homozygous
   */
  @Test
  void testCyp2c19_s1s2(TestInfo testInfo) throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("CYP2C19", "rs12769205", "A", "G")
        .variation("CYP2C19", "rs4244285", "G", "A")
        .variation("CYP2C19", "rs3758581", "G", "G");
    testWrapper.execute(s_otherOutsideCallFilePath);

    testWrapper.testCalledByMatcher("CYP2C19");
    testWrapper.testPrintCpicCalls( "CYP2C19", "*1/*2");

    testWrapper.testMatchedAnnotations("amitriptyline", 2);
    testWrapper.testAnyMatchFromSource("amitriptyline", DataSource.CPIC);
    testWrapper.testAnyMatchFromSource("amitriptyline", DataSource.DPWG);
    testWrapper.testMatchedAnnotations("citalopram", 2);
    testWrapper.testAnyMatchFromSource("citalopram", DataSource.CPIC);
    testWrapper.testAnyMatchFromSource("citalopram", DataSource.DPWG);
    testWrapper.testMatchedAnnotations("clomipramine", 3);
    testWrapper.testAnyMatchFromSource("clomipramine", DataSource.CPIC);
    testWrapper.testAnyMatchFromSource("clomipramine", DataSource.DPWG);
    testWrapper.testMatchedAnnotations("ivacaftor", 0);

    GeneReport cyp2c19report = testWrapper.getContext().getGeneReport(DataSource.CPIC, "CYP2C19");

    // make sure the variant in question is not a het call
    VariantReport vr = cyp2c19report.findVariantReport("rs58973490")
        .orElseThrow(() -> new RuntimeException("Variant missing from test data"));
    assertFalse(vr.isHetCall());

    // the variant is hom so ambiguity message should not apply and, thus, no matching messages
    assertEquals(0, cyp2c19report.getMessages().stream()
        .filter(m -> m.getExceptionType().equals(MessageAnnotation.TYPE_AMBIGUITY) &&
            Objects.requireNonNull(m.getMatches().getVariant()).equals("rs58973490"))
        .count());

    testWrapper.testAnyMatchFromSource("amitriptyline", DataSource.CPIC);
    // the variant is hom so ambiguity message should not match
    testWrapper.testMessageCountForDrug(DataSource.CPIC, "amitriptyline", 0);

    // CYP2C19 reference is *38, not *1, so should not have reference message
    testWrapper.testMessageCountForGene(DataSource.CPIC, "CYP2C19", 0);
  }

  @Test
  void testClomipramineCall(TestInfo testInfo) throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("CYP2C19", "rs12769205", "G", "G")
        .variation("CYP2C19", "rs4244285", "A", "A")
        .variation("CYP2C19", "rs3758581", "G", "G");
    testWrapper.execute(s_otherOutsideCallFilePath);

    testWrapper.testCalledByMatcher("CYP2C19");
    testWrapper.testPrintCpicCalls( "CYP2C19", "*2/*2");

    testWrapper.testMatchedAnnotations("amitriptyline", 2);
    testWrapper.testAnyMatchFromSource("amitriptyline", DataSource.CPIC);
    testWrapper.testAnyMatchFromSource("amitriptyline", DataSource.DPWG);
    testWrapper.testMatchedAnnotations("clomipramine", 3);
    testWrapper.testAnyMatchFromSource("clomipramine", DataSource.CPIC);
    testWrapper.testAnyMatchFromSource("clomipramine", DataSource.DPWG);
    testWrapper.testMatchedAnnotations("desipramine", 1);
    testWrapper.testAnyMatchFromSource("desipramine", DataSource.CPIC);
    testWrapper.testMatchedAnnotations("doxepin", 2);
    testWrapper.testAnyMatchFromSource("doxepin", DataSource.CPIC);
    testWrapper.testAnyMatchFromSource("doxepin", DataSource.DPWG);
    testWrapper.testMatchedAnnotations("imipramine", 3);
    testWrapper.testAnyMatchFromSource("imipramine", DataSource.CPIC);
    testWrapper.testAnyMatchFromSource("imipramine", DataSource.DPWG);
    testWrapper.testMatchedAnnotations("nortriptyline", 2);
    testWrapper.testAnyMatchFromSource("nortriptyline", DataSource.CPIC);
    testWrapper.testAnyMatchFromSource("nortriptyline", DataSource.DPWG);
    testWrapper.testMatchedAnnotations("trimipramine", 1);

    testWrapper.testMatchedAnnotations("clopidogrel", 4);
    testWrapper.testAnyMatchFromSource("clopidogrel", DataSource.CPIC);
    testWrapper.testAnyMatchFromSource("clopidogrel", DataSource.DPWG);

    testWrapper.testMatchedAnnotations("lansoprazole", 2);
    testWrapper.testAnyMatchFromSource("lansoprazole", DataSource.CPIC);
    testWrapper.testAnyMatchFromSource("lansoprazole", DataSource.DPWG);

    // voriconazole has 2 populations with recommendations so should have 2 matching annotations from CPIC
    // and 1 from DPWG
    testWrapper.testMatchedAnnotations("voriconazole", DataSource.CPIC, 2);
    testWrapper.testMatchedAnnotations("voriconazole", DataSource.DPWG, 1);
  }

  @Test
  void testCyp2c19noCall(TestInfo testInfo) throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("CYP2C19", "rs12769205", "A", "G")
        .variation("CYP2C19", "rs4244285", "A", "A");
    testWrapper.execute(s_otherOutsideCallFilePath);

    testWrapper.testNotCalledByMatcher("CYP2C19");

    testWrapper.testNoMatchFromSource("citalopram", DataSource.CPIC);
    testWrapper.testNoMatchFromSource("citalopram", DataSource.DPWG);
    testWrapper.testNoMatchFromSource("ivacaftor", DataSource.CPIC);
    testWrapper.testNoMatchFromSource("ivacaftor", DataSource.DPWG);
  }

  @Test
  void testCyp2c19s4bs17rs28399504missing(TestInfo testInfo) throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("CYP2C19", "rs12248560", "T", "T")
        .missing("CYP2C19", "rs28399504")
        .variation("CYP2C19", "rs3758581", "G", "G");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("CYP2C19");
    testWrapper.testPrintCpicCalls("CYP2C19", "*4/*4", "*4/*17", "*17/*17");

    testWrapper.testMatchedAnnotations("citalopram", 6);
    testWrapper.testMatchedAnnotations("citalopram", DataSource.CPIC, 3);
    testWrapper.testMatchedAnnotations("citalopram", DataSource.DPWG, 3);
    testWrapper.testAnyMatchFromSource("citalopram", DataSource.CPIC);
    testWrapper.testAnyMatchFromSource("citalopram", DataSource.DPWG);
  }

  @Test
  void testCyp2c19s1s4het(TestInfo testInfo) throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("CYP2C19", "rs12248560", "T", "T")
        .variation("CYP2C19", "rs28399504", "A", "G")
        .variation("CYP2C19", "rs3758581", "G", "G");
    testWrapper.execute(s_outsideCallFilePath);

    testWrapper.testCalledByMatcher("CYP2C19");
    testWrapper.testReportable("CYP2D6");

    testWrapper.testPrintCpicCalls("CYP2D6", "*1/*4");
    testWrapper.testPrintCpicCalls("CYP2C19", "*4/*17");

    assertTrue(testWrapper.getContext().getGeneReport(DataSource.CPIC, "CYP2D6").isOutsideCall());
  }

  @Test
  void testCyp2c19s1s4missingS1(TestInfo testInfo) throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("CYP2C19", "rs12248560", "C", "T")
        .variation("CYP2C19", "rs28399504", "A", "G")
        .missing("CYP2C19", "rs3758581");
    testWrapper.execute(s_outsideCallFilePath);

    testWrapper.testCalledByMatcher("CYP2C19");
    testWrapper.testReportable("CYP2D6");

    testWrapper.testPrintCpicCalls("CYP2D6", "*1/*4");
    testWrapper.testPrintCpicCalls("CYP2C19",  "*1/*4", "*4/*38");

    assertTrue(testWrapper.getContext().getGeneReport(DataSource.CPIC, "CYP2D6").isOutsideCall());
    GeneReport cyp2c19report = testWrapper.getContext().getGeneReport(DataSource.CPIC, "CYP2C19");
    assertTrue(cyp2c19report.isMissingVariants());

    assertFalse(cyp2c19report.isPhased());
    assertTrue(cyp2c19report.findVariantReport("rs12248560").map(VariantReport::isHetCall).orElse(false));
    assertTrue(cyp2c19report.findVariantReport("rs3758581").map(VariantReport::isMissing).orElse(false));

    assertTrue(cyp2c19report.hasHaplotype("*38"));

    // message is for *1/*4 being ambiuguous with unphased data
    assertEquals(2, cyp2c19report.getMessages().stream()
        .filter(m -> m.getExceptionType().equals(MessageAnnotation.TYPE_AMBIGUITY))
        .count());

    testWrapper.testAnyMatchFromSource("amitriptyline", DataSource.CPIC);
    testWrapper.testMessageCountForDrug(DataSource.CPIC, "amitriptyline", 2);
  }

  @Test
  void testCyp2c19s4s17(TestInfo testInfo) throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("CYP2C19", "rs12248560", "C", "T")
        .variation("CYP2C19", "rs28399504", "A", "G")
        .variation("CYP2C19", "rs3758581", "G", "G");
    testWrapper.execute(s_outsideCallFilePath);

    testWrapper.testCalledByMatcher("CYP2C19");
    testWrapper.testReportable("CYP2D6");

    testWrapper.testPrintCpicCalls("CYP2D6", "*1/*4");
    testWrapper.testPrintCpicCalls("CYP2C19", "*1/*4");

    assertTrue(testWrapper.getContext().getGeneReport(DataSource.CPIC, "CYP2D6").isOutsideCall());
  }

  @Test
  void testCftrRefRef(TestInfo testInfo) throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("CFTR");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("CFTR");
    testWrapper.testPrintCpicCalls("CFTR", "No CPIC variants found");
    testWrapper.testLookup("CFTR", "ivacaftor non-responsive CFTR sequence", "ivacaftor non-responsive CFTR sequence");

    testWrapper.testMatchedAnnotations("ivacaftor", 1);
  }

  @Test
  void testCftrD1270NHet(TestInfo testInfo) throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("CFTR", "rs11971167", "G", "A");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("CFTR");
    testWrapper.testPrintCpicCalls("CFTR", "D1270N (heterozygous)");
    testWrapper.testLookup("CFTR", "ivacaftor non-responsive CFTR sequence", "D1270N");

    testWrapper.testMatchedAnnotations("ivacaftor", 1);
  }

  @Test
  void testCftrD1270NG551D(TestInfo testInfo) throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("CFTR", "rs11971167", "G", "A")
        .variation("CFTR", "rs75527207", "G", "A");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("CFTR");
    testWrapper.testPrintCpicCalls("CFTR", "D1270N/G551D");
    testWrapper.testLookup("CFTR", "G551D", "D1270N");

    testWrapper.testMatchedAnnotations("ivacaftor", 1);
  }

  @Test
  void testRosuvastatin(TestInfo testInfo) throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("ABCG2", "rs2231142", "G", "T")
        .variation("SLCO1B1", "rs56101265", "T", "C");
    Path vcfFile = testWrapper.execute(null);

    testWrapper.testCalledByMatcher("ABCG2", "SLCO1B1");
    testWrapper.testPrintCpicCalls("SLCO1B1", "*1/*2");

    testWrapper.testMatchedAnnotations("rosuvastatin", 1);

    // no dpyd - should not have DPYD warning
    Path reporterOutput = vcfFile.getParent().resolve(PharmCAT.getBaseFilename(vcfFile) + ".report.html");
    Document document = Jsoup.parse(reporterOutput.toFile());
    Elements capecitabineSection = document.getElementsByClass("capecitabine");
    assertEquals(1, capecitabineSection.size());
    assertEquals(0, capecitabineSection.get(0).getElementsByClass("alert-info").size());

    Elements dpydSection = document.select(".gene.dpyd");
    assertEquals(1, dpydSection.size());
    assertEquals(1, dpydSection.get(0).getElementsByClass("no-data").size());
  }

  @Test
  void testAmitryptylineCallWoCyp2c19(TestInfo testInfo) throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("DPYD");
    testWrapper.execute(s_outsideCallFilePath);

    testWrapper.testReportable("CYP2D6");
    testWrapper.testPrintCpicCalls("CYP2D6", "*1/*4");
    testWrapper.testLookup("CYP2D6", "*1", "*4");

    testWrapper.testMatchedAnnotations("amitriptyline", 2);
    testWrapper.testAnyMatchFromSource("amitriptyline", DataSource.CPIC);
    testWrapper.testAnyMatchFromSource("amitriptyline", DataSource.DPWG);
  }

  @Test
  void testSlco1b1HomWild(TestInfo testInfo) throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("SLCO1B1");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("SLCO1B1");
    testWrapper.testPrintCpicCalls("SLCO1B1", "*1/*1");
    testWrapper.testLookup("SLCO1B1", "*1", "*1");

    DrugReport drugReport = testWrapper.getContext().getDrugReport(DataSource.CPIC, "simvastatin");
    assertNotNull(drugReport);
    GuidelineReport guidelineReport = drugReport.getGuidelines().get(0);
    assertEquals(1, guidelineReport.getAnnotations().size());
    AnnotationReport annotationReport = guidelineReport.getAnnotations().get(0);
    assertTrue(annotationReport.getHighlightedVariants().contains("rs4149056:T/T"));
  }

  @Test
  void testSlco1b1HomVar(TestInfo testInfo) throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("SLCO1B1", "rs2306283", "A", "G")
        .variation("SLCO1B1", "rs4149056", "C", "C");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("SLCO1B1");
    testWrapper.testPrintCpicCalls("SLCO1B1", "*5/*15");
    testWrapper.testLookup("SLCO1B1", "*5", "*15");
  }

  @Test
  void testSlco1b1Test5(TestInfo testInfo) throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("SLCO1B1", "rs2306283", "A", "G")
        .variation("SLCO1B1", "rs11045852", "A", "G")
        .variation("SLCO1B1", "rs74064213", "A", "G");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("SLCO1B1");
    testWrapper.testPrintCpicCalls("SLCO1B1", "*1/*44");
    testWrapper.testLookup("SLCO1B1", "*1", "*44");
  }

  @Test
  void testSlco1b1Test3(TestInfo testInfo) throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("SLCO1B1", "rs2306283", "A", "G")
        .variation("SLCO1B1", "rs4149056", "T", "C");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("SLCO1B1");
    testWrapper.testPrintCpicCalls("SLCO1B1", "*1/*15");
    testWrapper.testLookup("SLCO1B1", "*1", "*15");
  }

  @Test
void testSlco1b1Test4(TestInfo testInfo) throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("SLCO1B1", "rs4149056", "T", "C")
        .variation("SLCO1B1", "rs71581941", "C", "T");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("SLCO1B1");
    testWrapper.testPrintCpicCalls("SLCO1B1", "*5/*45");
    testWrapper.testLookup("SLCO1B1", "*5", "*45");

    testWrapper.testMatchedAnnotations("simvastatin", 1);
  }

  @Test
  void testDpydS1S2B(TestInfo testInfo) throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .phased()
        .variation("DPYD", "rs3918290", "C", "T")
        .variation("DPYD", "rs1801159", "C", "T");
    Path vcfFile = testWrapper.execute(null);

    testWrapper.testCalledByMatcher("DPYD");

    GeneReport dpyd = testWrapper.getContext().getGeneReport(DataSource.CPIC, "DPYD");
    testWrapper.testPrintCpicCalls("DPYD", "c.1627A>G (*5)/c.1905+1G>A (*2A)");
    testWrapper.testLookup("DPYD", "c.1627A>G (*5)", "c.1905+1G>A (*2A)");

    testWrapper.testMatchedAnnotations("fluorouracil", 2);
    testWrapper.testMatchedAnnotations("capecitabine", 2);

    // should have DPYD warning
    Path reporterOutput = vcfFile.getParent().resolve(PharmCAT.getBaseFilename(vcfFile) + ".report.html");
    Document document = Jsoup.parse(reporterOutput.toFile());
    Elements capecitabineSection = document.getElementsByClass("capecitabine");
    assertEquals(1, capecitabineSection.size());
    Elements capecitabineMsgs = capecitabineSection.get(0).getElementsByClass("alert-info");
    assertEquals(1, capecitabineMsgs.size());
    assertTrue(capecitabineMsgs.get(0).text().contains("lowest activity"));

    Elements dpydSection = document.select(".gene.dpyd");
    assertEquals(1, dpydSection.size());
    assertEquals(0, dpydSection.get(0).getElementsByClass("no-data").size());
  }

  @Test
  void testDpydUnphased(TestInfo testInfo) throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("DPYD", "rs3918290", "C", "T")
        .variation("DPYD", "rs1801159", "C", "T");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("DPYD");

    GeneReport dpyd = testWrapper.getContext().getGeneReport(DataSource.CPIC, "DPYD");
    assertEquals(1, dpyd.getRecommendationDiplotypes().size());
    testWrapper.testPrintCpicCalls("DPYD", "c.1627A>G (*5)", "c.1905+1G>A (*2A)");
    testWrapper.testLookupByActivity("DPYD", "1.0");

    testWrapper.testMatchedAnnotations("fluorouracil", 2);
    testWrapper.testMatchedAnnotations("capecitabine", 2);
  }

  @Test
  void testDpydUnphasedMultiple(TestInfo testInfo) throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("DPYD", "rs183385770", "C", "T")  // 0 activity value
        .variation("DPYD", "rs186169810", "A", "C") // 0.5 activity value
        .variation("DPYD", "rs112766203", "G", "A") // 0.5 activity value
        .variation("DPYD", "rs144395748", "G", "C"); // 1 activity value
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("DPYD");

    GeneReport dpyd = testWrapper.getContext().getGeneReport(DataSource.CPIC, "DPYD");
    assertEquals(1, dpyd.getRecommendationDiplotypes().size());
    testWrapper.testPrintCpicCalls("DPYD", "c.1024G>A", "c.1314T>G", "c.1358C>G", "c.2279C>T");
    testWrapper.testLookupByActivity("DPYD", "0.5");

    testWrapper.testMatchedAnnotations("fluorouracil", 1);
    testWrapper.testMatchedAnnotations("capecitabine", 1);
  }

  @Test
  void testDpydC2846het(TestInfo testInfo) throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .phased()
        .variation("DPYD", "rs67376798", "T", "A");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("DPYD");
    testWrapper.testPrintCpicCalls("DPYD", "Reference/c.2846A>T");
    testWrapper.testLookup("DPYD", "Reference", "c.2846A>T");

    testWrapper.testAnyMatchFromSource("fluorouracil", DataSource.CPIC);
    testWrapper.testAnyMatchFromSource("fluorouracil", DataSource.DPWG);

    testWrapper.testAnyMatchFromSource("capecitabine", DataSource.CPIC);
    testWrapper.testAnyMatchFromSource("capecitabine", DataSource.DPWG);
  }

  /**
   * This test puts 2 alleles on each strand of a phased DPYD and then asserts that the least-function allele is used
   * for lookup on each of the strands.
   */
  @Test
  void testDpydPhasedMultiTrans(TestInfo testInfo) throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .phased()
        .variation("DPYD", "rs67376798", "A", "T") // Strand 1 decreased - c.2846A>T
        .variation("DPYD", "rs72547601", "C", "T") // Strand 1 no function - c.2933A>G
        .variation("DPYD", "rs60139309", "T", "C") // Strand 2 normal function - c.2582A>G
        .variation("DPYD", "rs139834141", "C", "T") // Strand 2 normal function - c.498G>A
    ;
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("DPYD");
    testWrapper.testMatcher("DPYD", "c.498G>A + c.2582A>G/c.2846A>T + c.2933A>G");
    testWrapper.testPrintCpicCalls("DPYD", "c.498G>A + c.2582A>G/c.2846A>T + c.2933A>G");
    assertEquals(1, testWrapper.getContext().getGeneReport(DataSource.CPIC, "DPYD").getRecommendationDiplotypes().size());
    testWrapper.testLookup("DPYD", "c.498G>A", "c.2933A>G");

    testWrapper.testAnyMatchFromSource("fluorouracil", DataSource.CPIC);

    testWrapper.testAnyMatchFromSource("capecitabine", DataSource.CPIC);
  }

  /**
   * This test is the same as the previous test but DPYD is unphased instead of phased. This means the individual found
   * alleles should be reported and then the two least-function alleles should be used for recommendation lookup.
   */
  @Test
  void testDpydUnphasedMulti(TestInfo testInfo) throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("DPYD", "rs67376798", "A", "T") // decreased - c.2846A>T
        .variation("DPYD", "rs72547601", "C", "T") // no function - c.2933A>G
        .variation("DPYD", "rs60139309", "T", "C") // normal function - c.2582A>G
        .variation("DPYD", "rs139834141", "C", "T") // normal function - c.498G>A
    ;
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("DPYD");
    testWrapper.testPrintCpicCalls("DPYD", "c.498G>A", "c.2582A>G", "c.2846A>T", "c.2933A>G");
    testWrapper.testLookup("DPYD", "c.2933A>G", "c.2846A>T");

    testWrapper.testAnyMatchFromSource("fluorouracil", DataSource.CPIC);

    testWrapper.testAnyMatchFromSource("capecitabine", DataSource.CPIC);
  }

  @Test
  void testDpydS12het(TestInfo testInfo) throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .phased()
        .variation("DPYD", "rs78060119", "C", "A");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("DPYD");
    testWrapper.testPrintCpicCalls("DPYD", "Reference/c.1156G>T (*12)");
    testWrapper.testLookup("DPYD", "Reference", "c.1156G>T (*12)");

    testWrapper.testMatchedAnnotations("fluorouracil", 1);
    testWrapper.testMatchedAnnotations("capecitabine", 1);
  }

  @Test
  void testDpydHomNoFunctionEffectivelyPhased(TestInfo testInfo) throws Exception {
    // effectively phased
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("DPYD")
        .variation("DPYD", "rs72549310", "A", "A")   // c.61C>T, hom variant (No function)
        .variation("DPYD", "rs150385342", "C", "T"); // c.313G>A het variant (Normal function)
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("DPYD");
    testWrapper.testMatcher("DPYD", "c.61C>T/c.61C>T + c.313G>A");
    testWrapper.testPrintCpicCalls("DPYD", "c.61C>T/c.61C>T + c.313G>A");
    assertEquals(1, testWrapper.getContext().getGeneReport(DataSource.CPIC, "DPYD").getRecommendationDiplotypes().size());
    testWrapper.testLookup("DPYD", "c.61C>T", "c.61C>T");

    testWrapper.testMatchedAnnotations("fluorouracil", 1);
    testWrapper.testMatchedAnnotations("capecitabine", 1);
  }

  @Test
  void testDpydHomNoFunctionPhased(TestInfo testInfo) throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .phased()
        .reference("DPYD")
        .variation("DPYD", "rs72549310", "A", "A")   // c.61C>T, hom variant (No function)
        .variation("DPYD", "rs150385342", "C", "T"); // c.313G>A het variant (Normal function)
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("DPYD");
    testWrapper.testMatcher("DPYD", "c.61C>T/c.61C>T + c.313G>A");
    testWrapper.testPrintCpicCalls("DPYD", "c.61C>T/c.61C>T + c.313G>A");
    assertEquals(1, testWrapper.getContext().getGeneReport(DataSource.CPIC, "DPYD").getRecommendationDiplotypes().size());
    testWrapper.testLookup("DPYD", "c.61C>T", "c.61C>T");

    testWrapper.testMatchedAnnotations("fluorouracil", 1);
    testWrapper.testMatchedAnnotations("capecitabine", 1);
  }

  @Test
  void testDpydHomNoFunctionUnphased(TestInfo testInfo) throws Exception {
    // effectively phased
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("DPYD")
        .variation("DPYD", "rs72547601", "C", "C") // c.2933A>G - no function
        .variation("DPYD", "rs67376798", "A", "T") // c.2846A>T - decreased
        .variation("DPYD", "rs60139309", "T", "C") // c.2582A>G - normal
    ;
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("DPYD");
    testWrapper.testMatcher("DPYD", "c.2582A>G", "c.2846A>T", "c.2933A>G");
    testWrapper.testPrintCpicCalls("DPYD", "c.2582A>G", "c.2846A>T", "c.2933A>G");
    assertEquals(1, testWrapper.getContext().getGeneReport(DataSource.CPIC, "DPYD").getRecommendationDiplotypes().size());
    testWrapper.testLookup("DPYD", "c.2933A>G", "c.2933A>G");

    testWrapper.testMatchedAnnotations("fluorouracil", 1);
    testWrapper.testMatchedAnnotations("capecitabine", 1);
  }


  /**
   * This tests a special case of SLCO1B1. The gene in this scenario is "uncalled" by the matcher due the sample VCF
   * data. However, SLCO1B1 has an override that will display the rs4149056 diplotype regardless of call status. That
   * same override will assign alleles to use for recommendation lookup
   */
  @Test
  void testSlco1b1TestMulti(TestInfo testInfo) throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("SLCO1B1", "rs2306283", "G", "G")
        .variation("SLCO1B1", "rs4149056", "T", "C")
        .variation("SLCO1B1", "rs11045853", "A", "A")
        .variation("SLCO1B1", "rs72559748", "G", "G");
    testWrapper.execute(null);

    testWrapper.testNotCalledByMatcher("SLCO1B1");
    testWrapper.testInferredCalls("SLCO1B1", "rs4149056T/rs4149056C");
    testWrapper.testLookup("SLCO1B1", "*1", "*5");

    testWrapper.testMatchedAnnotations("simvastatin", 1);
  }

  @Test
  void testUgt1a1Phased(TestInfo testInfo) throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .phased()
        .variation("UGT1A1", "rs887829", "C", "T");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("UGT1A1");
    testWrapper.testPrintCpicCalls("UGT1A1", "*1/*80");
    testWrapper.testLookup("UGT1A1", "*1", "*80");
  }

  @Test
  void testUgt1a1Unphased(TestInfo testInfo) throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .phased()
        .variation("UGT1A1", "rs887829", "C", "T");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("UGT1A1");
    testWrapper.testPrintCpicCalls("UGT1A1", "*1/*80");
    testWrapper.testLookup("UGT1A1", "*1", "*80");
  }

  @Test
  void testUgt1a1s1s1(TestInfo testInfo) throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("UGT1A1");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("UGT1A1");
    testWrapper.testPrintCpicCalls("UGT1A1", "*1/*1");
    testWrapper.testLookup("UGT1A1", "*1", "*1");
  }

  @Test
  void testUgt1a1S1S80S28(TestInfo testInfo) throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("UGT1A1", "rs887829", "C", "T")
        .variation("UGT1A1", "rs3064744", "TA(7)", "TA(8)");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("UGT1A1");
    testWrapper.testPrintCpicCalls("UGT1A1", "*1/*80+*28");
    testWrapper.testLookup("UGT1A1", "*1", "*80+*28");
  }

  @Test
  void testUgt1a1S28S37(TestInfo testInfo) throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("UGT1A1", "rs3064744", "TA(8)", "TA(9)");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("UGT1A1");
    testWrapper.testPrintCpicCalls("UGT1A1", "*28/*37");
    testWrapper.testLookup("UGT1A1", "*37", "*28");
  }

  @Test
  void testUgt1a1s28s80phased(TestInfo testInfo) throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .phased()
        .variation("UGT1A1", "rs887829", "C", "T")
        .variation("UGT1A1", "rs3064744", "TA(7)", "TA(8)");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("UGT1A1");
    testWrapper.testPrintCpicCalls("UGT1A1", "*1/*80+*28");
    testWrapper.testLookup("UGT1A1", "*1", "*80+*28");

    // the guideline should not have an ambiguity message
    testWrapper.testMatchedAnnotations("atazanavir", 1);
    testWrapper.testMessageCountForDrug(DataSource.CPIC, "atazanavir", 0);

    testWrapper.testMessageCountForGene(DataSource.CPIC, "UGT1A1", 2);
    testWrapper.testGeneHasMessage(DataSource.CPIC, "UGT1A1", "reference-allele");
  }

  @Test
  void testUgt1a1s28s80s6s60phased(TestInfo testInfo) throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .phased()
        .variation("UGT1A1", "rs887829", "T", "C")
        .variation("UGT1A1", "rs3064744", "TA(8)", "TA(7)")
        .variation("UGT1A1", "rs4148323", "G", "A");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("UGT1A1");
    testWrapper.testPrintCpicCalls("UGT1A1", "*6/*80+*28");
    testWrapper.testLookup("UGT1A1", "*6", "*80+*28");
  }

  @Test
  void testUgt1a1s28s80s6s60unphased(TestInfo testInfo) throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("UGT1A1", "rs887829", "T", "C")
        .variation("UGT1A1", "rs3064744", "TA(8)", "TA(7)")
        .variation("UGT1A1", "rs4148323", "G", "A");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("UGT1A1");
    testWrapper.testPrintCpicCalls("UGT1A1", "*6/*80+*28");
    testWrapper.testLookup("UGT1A1", "*6", "*80+*28");
  }

  @Test
  void testUgt1a1s6s6(TestInfo testInfo) throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("UGT1A1", "rs4148323", "A", "A");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("UGT1A1");
    testWrapper.testPrintCpicCalls("UGT1A1", "*6/*6");
    testWrapper.testLookup("UGT1A1", "*6", "*6");
  }

  @Test
  void testUgt1a1s6s60s80s28MissingPhased(TestInfo testInfo) throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .phased()
        .missing("UGT1A1", "rs3064744")
        .variation("UGT1A1", "rs887829", "C", "T")
        .variation("UGT1A1", "rs4148323", "A", "G");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("UGT1A1");
    testWrapper.testPrintCpicCalls("UGT1A1", "*6/*80", "*6/*80+*28", "*6/*80+*37");
    testWrapper.testLookup("UGT1A1", "*6", "*80");
    testWrapper.testLookup("UGT1A1", "*6", "*80+*28");
    testWrapper.testLookup("UGT1A1", "*6", "*80+*37");
  }

  @Test
  void testUgt1a1s6s60s80s28MissingUnphased(TestInfo testInfo) throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .missing("UGT1A1", "rs3064744")
        .variation("UGT1A1", "rs887829", "C", "T")
        .variation("UGT1A1", "rs4148323", "A", "G");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("UGT1A1");
    testWrapper.testPrintCpicCalls("UGT1A1", "*6/*80", "*6/*80+*28", "*6/*80+*37");
    testWrapper.testLookup("UGT1A1", "*6", "*80");
    testWrapper.testLookup("UGT1A1", "*6", "*80+*28");
    testWrapper.testLookup("UGT1A1", "*6", "*80+*37");
  }

  @Test
  void testUgt1a1s80s28missing(TestInfo testInfo) throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .missing("UGT1A1", "rs3064744")
        .variation("UGT1A1", "rs887829", "C", "T");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("UGT1A1");
    testWrapper.testPrintCpicCalls("UGT1A1", "*1/*80", "*1/*80+*28", "*1/*80+*37");
    testWrapper.testLookup("UGT1A1", "*1", "*80");
    testWrapper.testLookup("UGT1A1", "*1", "*80+*28");
    testWrapper.testLookup("UGT1A1", "*1", "*80+*37");
  }

  @Test
  void testUgt1a1na12717(TestInfo testInfo) throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("UGT1A1", "rs887829", "T", "T")
        .variation("UGT1A1", "rs3064744", "TA(8)", "TA(7)");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("UGT1A1");
    testWrapper.testPrintCpicCalls("UGT1A1", "*80/*80+*28");
    testWrapper.testLookup("UGT1A1", "*80", "*80+*28");
  }

  @Test
  void testUgt1a1s28homMissing(TestInfo testInfo) throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .missing("UGT1A1", "rs887829")
        .variation("UGT1A1", "rs3064744", "TA(8)", "TA(8)");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("UGT1A1");
    testWrapper.testPrintCpicCalls("UGT1A1", "*28/*28", "*28/*80+*28", "*80+*28/*80+*28");
    testWrapper.testLookup("UGT1A1", "*28", "*28");
    testWrapper.testLookup("UGT1A1", "*28", "*80+*28");
    testWrapper.testLookup("UGT1A1", "*80+*28", "*80+*28");
  }

  @Test
  void testUgt1a1s28s60Hom(TestInfo testInfo) throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("UGT1A1", "rs3064744", "TA(8)", "TA(7)");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("UGT1A1");
    testWrapper.testPrintCpicCalls("UGT1A1", "*1/*28");
    testWrapper.testLookup("UGT1A1", "*1", "*28");
  }

  @Test
  void testUgt1a1s27s28unphaseds80s60missing(TestInfo testInfo) throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .missing("UGT1A1", "rs887829")
        .variation("UGT1A1", "rs3064744", "TA(8)", "TA(7)")
        .variation("UGT1A1", "rs35350960", "C", "A");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("UGT1A1");
    testWrapper.testPrintCpicCalls("UGT1A1", "*27/*28", "*27/*80+*28");
    testWrapper.testLookup("UGT1A1", "*27", "*28");
    testWrapper.testLookup("UGT1A1", "*27", "*80+*28");
  }

  @Test
  void testUgt1a1HG00436(TestInfo testInfo) throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .phased()
        .variation("UGT1A1", "rs887829", "T", "C")
        .variation("UGT1A1", "rs3064744", "TA(8)", "TA(7)")
        .variation("UGT1A1", "rs35350960", "A", "C");
    testWrapper.execute(null);

    testWrapper.testNotCalledByMatcher("UGT1A1");
  }

  @Test
  void testUgt1a1s1s80s27s60s28missingphased(TestInfo testInfo) throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .phased()
        .missing("UGT1A1", "rs3064744")
        .variation("UGT1A1", "rs887829", "T", "C")
        .variation("UGT1A1", "rs35350960", "A", "C");
    testWrapper.execute(null);

    testWrapper.testNotCalledByMatcher("UGT1A1");
    GeneReport geneReport = testWrapper.getContext().getGeneReport(DataSource.CPIC, "UGT1A1");
    assertTrue(geneReport.isPhased());
  }

  @Test
  void testUgt1a1s1s60s80s6phased(TestInfo testInfo) throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .phased()
        .variation("UGT1A1", "rs887829", "T", "C")
        .variation("UGT1A1", "rs35350960", "A", "C");
    testWrapper.execute(null);

    testWrapper.testNotCalledByMatcher("UGT1A1");
    GeneReport geneReport = testWrapper.getContext().getGeneReport(DataSource.CPIC, "UGT1A1");
    assertTrue(geneReport.isPhased());
  }

  @Test
  void testUgt1a1s1s60s80s28s6phased(TestInfo testInfo) throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .phased()
        .variation("UGT1A1", "rs887829", "T", "C")
        .variation("UGT1A1", "rs3064744", "TA(8)", "TA(7)")
        .variation("UGT1A1", "rs35350960", "A", "C");
    testWrapper.execute(null);

    testWrapper.testNotCalledByMatcher("UGT1A1");
    GeneReport geneReport = testWrapper.getContext().getGeneReport(DataSource.CPIC, "UGT1A1");
    assertTrue(geneReport.isPhased());
  }

  @Test
  void testUgt1a1s1s37s80s60phased(TestInfo testInfo) throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .phased()
        .variation("UGT1A1", "rs887829", "T", "C")
        .variation("UGT1A1", "rs3064744", "TA(9)", "TA(7)");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("UGT1A1");
    testWrapper.testReportable("UGT1A1");
    testWrapper.testPrintCpicCalls("UGT1A1", "*1/*80+*37");
    testWrapper.testLookup("UGT1A1", "*1", "*80+*37");
    GeneReport geneReport = testWrapper.getContext().getGeneReport(DataSource.CPIC, "UGT1A1");
    assertTrue(geneReport.isPhased());
  }

  @Test
  void testCyp3a5Missing3Message(TestInfo testInfo) throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .missing("CYP3A5", "rs776746");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("CYP3A5");
    testWrapper.testReportable("CYP3A5");
    testWrapper.testPrintCpicCalls("CYP3A5", "*1/*1");
    testWrapper.testLookup("CYP3A5", "*1", "*1");

    GeneReport gene = testWrapper.getContext().getGeneReport(DataSource.CPIC, "CYP3A5");
    // rs776746 should be missing from this report
    assertNotNull(gene.getVariantReports());
    assertTrue(gene.getVariantReports().stream().anyMatch(v -> v.isMissing() && v.getDbSnpId().equals("rs776746")));

    // the guideline should have a matching message
    assertTrue(testWrapper.getContext().getDrugReports().get(DataSource.CPIC).values().stream()
        .filter(r -> r.getName().equals("tacrolimus"))
        .allMatch(r -> r.getMessages().size() > 0));

    assertFalse(gene.isPhased());
  }

  @Test
  void testCyp3a5v1(TestInfo testInfo) throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("CYP3A5", "rs776746", "T", "C");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("CYP3A5");
    testWrapper.testReportable("CYP3A5");
    testWrapper.testPrintCpicCalls("CYP3A5", "*1/*3");
    testWrapper.testLookup("CYP3A5", "*1", "*3");
  }

  @Test
  void testCyp3a5v2(TestInfo testInfo) throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("CYP3A5", "rs28383479", "C", "T")
        .variation("CYP3A5", "rs776746", "C", "T")
    ;
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("CYP3A5");
    testWrapper.testPrintCpicCalls("CYP3A5", "*3/*9");
    testWrapper.testLookup("CYP3A5", "*3", "*9");
  }

  @Test
  void testCyp3a5v3(TestInfo testInfo) throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("CYP3A5", "rs776746", "C", "C")
    ;
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("CYP3A5");
    testWrapper.testPrintCpicCalls("CYP3A5", "*3/*3");
    testWrapper.testLookup("CYP3A5", "*3", "*3");
  }

  @Test
  void testCyp3a5v4(TestInfo testInfo) throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("CYP3A5", "rs776746", "T", "C")
    ;
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("CYP3A5");
    testWrapper.testPrintCpicCalls("CYP3A5", "*1/*3");
    testWrapper.testLookup("CYP3A5", "*1", "*3");
  }

  @Test
  void testCyp3a5v5(TestInfo testInfo) throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("CYP3A5", "rs28383479", "T", "C")
        .variation("CYP3A5", "rs776746", "T", "C")
    ;
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("CYP3A5");
    testWrapper.testPrintCpicCalls("CYP3A5", "*3/*9");
    testWrapper.testLookup("CYP3A5", "*3", "*9");
  }

  @Test
  void testHlab(TestInfo testInfo) throws Exception {
    Path outsideCallPath = TestUtils.createTestFile(testInfo, ".tsv");
    try (PrintWriter writer = new PrintWriter(Files.newBufferedWriter(outsideCallPath))) {
      writer.println("HLA-B\t*15:02/*57:01");
    }

    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("CYP2C9");
    testWrapper.execute(outsideCallPath);

    testWrapper.testCalledByMatcher("CYP2C9");
    testWrapper.testNotCalledByMatcher("HLA-B");
    testWrapper.testReportable("CYP2C9");
    testWrapper.testReportable("HLA-B");
    testWrapper.testMatchedAnnotations("abacavir", DataSource.CPIC, 1);
    testWrapper.testMatchedAnnotations("abacavir", DataSource.DPWG, 1);
    testWrapper.testMatchedAnnotations("allopurinol", DataSource.CPIC, 1);
    testWrapper.testMatchedAnnotations("allopurinol", DataSource.DPWG, 1);
    testWrapper.testMatchedAnnotations("phenytoin", 4);
    testWrapper.testAnyMatchFromSource("phenytoin", DataSource.CPIC);
    testWrapper.testAnyMatchFromSource("phenytoin", DataSource.DPWG);
  }

  @Test
  void testHlabPhenotype(TestInfo testInfo) throws Exception {
    Path outsideCallPath = TestUtils.createTestFile(testInfo, ".tsv");
    try (PrintWriter writer = new PrintWriter(Files.newBufferedWriter(outsideCallPath))) {
      writer.println("HLA-B\t\t*57:01 positive");
    }

    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("CYP2C9");
    testWrapper.execute(outsideCallPath);

    testWrapper.testCalledByMatcher("CYP2C9");
    testWrapper.testNotCalledByMatcher("HLA-B");
    testWrapper.testReportable("CYP2C9");
    testWrapper.testReportable("HLA-B");
    testWrapper.testMatchedAnnotations("abacavir", 2);
    testWrapper.testMatchedAnnotations("abacavir", DataSource.CPIC, 1);
    testWrapper.testMatchedAnnotations("abacavir", DataSource.DPWG, 1);
    // allopurinol relies on a different allele for recs so no matches
    testWrapper.testMatchedAnnotations("allopurinol", 0);
    // phenytoin also relies on a different allele but there will be a match for DPWG since the recommendations are
    // split between the two genes on that side
    testWrapper.testMatchedAnnotations("phenytoin", 1);
    testWrapper.testNoMatchFromSource("phenytoin", DataSource.CPIC);
    testWrapper.testAnyMatchFromSource("phenytoin", DataSource.DPWG);
  }

  /**
   * An example report that shows a few different types of recommendation scenarios all in one report. The examples
   * shown are:
   * <ul>
   *   <li>celecoxib = 1 CPIC recommenation</li>
   *   <li>citalopram = 2 recommendations: 1 CPIC, 1 DPWG, 1 gene and it's called</li>
   *   <li>clomipramine = 2 recommendations: 1 CPIC, 1 DPWG, 2 gene but only 1 called</li>
   *   <li>carbamezepine = 3 CPIC recommendations on different populations</li>
   *   <li>clopidogrel = 4 recommendations: 3 CPIC on different pops, 1 DPWG</li>
   *   <li>flucloxacillin = 0 recommendations but the gene is reportable</li>
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

    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("CYP2C9")
        .variation("CYP2C19", "rs12769205", "G", "G")
        .variation("CYP2C19", "rs4244285", "A", "A")
        .variation("CYP2C19", "rs3758581", "G", "G");
    testWrapper.execute(outsideCallPath);

    testWrapper.testLookup("CYP2C19", "*2", "*2");
    testWrapper.testPrintCpicCalls("CYP2C19", "*2/*2");

    GeneReport cyp2c9 = testWrapper.getContext().getGeneReport(DataSource.CPIC, "CYP2C9");
    assertEquals(1, cyp2c9.getRecommendationDiplotypes().size());
    assertTrue(cyp2c9.getRecommendationDiplotypes().stream().allMatch(d -> d.getActivityScore().equals("2.0")));

    testWrapper.testReportable("CYP2C19", "CYP2C9", "HLA-A", "HLA-B");
    testWrapper.testMatchedAnnotations("celecoxib", 1);
    testWrapper.testAnyMatchFromSource("celecoxib", DataSource.CPIC);
    testWrapper.testMatchedAnnotations("citalopram", 2);
    testWrapper.testMatchedAnnotations("clomipramine", 2);
    testWrapper.testMatchedAnnotations("clopidogrel", 4);
    testWrapper.testMatchedAnnotations("clopidogrel", DataSource.CPIC, 3);
    testWrapper.testMatchedAnnotations("clopidogrel", DataSource.DPWG, 1);
    testWrapper.testNoMatchFromSource("flucloxacillin", DataSource.CPIC);
    testWrapper.testMatchedAnnotations("flucloxacillin", DataSource.DPWG, 1);
    testWrapper.testNoMatchFromSource("fluvoxamine", DataSource.CPIC);
    testWrapper.testNoMatchFromSource("fluvoxamine", DataSource.DPWG);
    testWrapper.testMatchedAnnotations("siponimod", 1);
    testWrapper.testAnyMatchFromSource("siponimod", DataSource.DPWG);

    testWrapper.testMatchedAnnotations("carbamazepine", 5);
  }

  @Test
  void testTpmtStar1s(TestInfo testInfo) throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("TPMT", "rs1800460", "C", "T")
        .variation("TPMT", "rs1142345", "T", "C");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("TPMT");
    testWrapper.testPrintCpicCalls("TPMT", "*1/*3A");
    testWrapper.testLookup("TPMT", "*1", "*3A");

    GeneReport tpmtReport = testWrapper.getContext().getGeneReport(DataSource.CPIC, "TPMT");
    assertEquals(43, tpmtReport.getVariantReports().size());
  }


  @Test
  void testCyp2c9star61(TestInfo testInfo) throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("CYP2C9", "rs1799853", "C", "T")
        .variation("CYP2C9", "rs202201137", "A", "G");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("CYP2C9");
    testWrapper.testPrintCpicCalls("CYP2C9", "*1/*61");
    testWrapper.testLookup("CYP2C9", "*1", "*61");
  }

  @Test
  void testCyp2c9star1Hom(TestInfo testInfo) throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("CYP2C9");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("CYP2C9");
    testWrapper.testPrintCpicCalls("CYP2C9", "*1/*1");
    testWrapper.testLookup("CYP2C9", "*1", "*1");
    testWrapper.testMatchedAnnotations("celecoxib", 1);
    testWrapper.testMatchedAnnotations("ibuprofen", 1);
    testWrapper.testMatchedAnnotations("lornoxicam", 1);
  }


  /**
   * Test CYP2B6 for a het *34 sample file. When doing the "top match" scenario this will only match to a 1/34 and,
   * thus, only match to a single recommendation. This test will have a different outcome when run in "all matches" mode
   * and should be compared with {@link #testCyp2b6star1star34AllMatch(TestInfo)}.
   */
  @Test
  void testCyp2b6star1star34(TestInfo testInfo) throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("CYP2B6", "rs34223104", "T", "C")
        .variation("CYP2B6", "rs3211371", "C", "A")
        .variation("CYP2B6", "rs3745274", "G", "T")
        .variation("CYP2B6", "rs2279343", "A", "G")
    ;
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("CYP2B6");
    testWrapper.testPrintCpicCalls("CYP2B6", "*1/*34");
    testWrapper.testLookup("CYP2B6", "*1", "*34");
    testWrapper.testMatchedAnnotations("efavirenz", 1);
  }

  /**
   * This test is just like {@link #testCyp2b6star1star34(TestInfo)} but run in "all matches" mode. This should result in 2
   * possible different calls coming from the matcher. These two have different phenotypes and, thus, match to different
   * recommendations.
   */
  @Test
  void testCyp2b6star1star34AllMatch(TestInfo testInfo) throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper(testInfo, true);
    testWrapper.getVcfBuilder()
        .variation("CYP2B6", "rs34223104", "T", "C")
        .variation("CYP2B6", "rs3211371", "C", "A")
        .variation("CYP2B6", "rs3745274", "G", "T")
        .variation("CYP2B6", "rs2279343", "A", "G")
    ;
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("CYP2B6");
    testWrapper.testPrintCpicCalls("CYP2B6", "*1/*34", "*33/*36");
    testWrapper.testLookup("CYP2B6", "*1", "*34");
    testWrapper.testLookup("CYP2B6", "*33", "*36");
    testWrapper.testMatchedAnnotations("efavirenz", 2);
  }


  /* NUDT15 */
  @Test
  void testNudt15Ref(TestInfo testInfo) throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("NUDT15");
    testWrapper.execute(null);

    testWrapper.testPrintCpicCalls("NUDT15", "*1/*1");
    testWrapper.testLookup("NUDT15", "*1", "*1");

    testWrapper.testMatchedAnnotations("azathioprine", 2);
    testWrapper.testMatchedAnnotations("mercaptopurine", 2);
    testWrapper.testAnyMatchFromSource("mercaptopurine", DataSource.CPIC);
    testWrapper.testAnyMatchFromSource("mercaptopurine", DataSource.DPWG);
    testWrapper.testMatchedAnnotations("thioguanine", 2);
    testWrapper.testAnyMatchFromSource("thioguanine", DataSource.CPIC);
    testWrapper.testAnyMatchFromSource("thioguanine", DataSource.DPWG);
  }

  @Test
  void testNudt15S2(TestInfo testInfo) throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("NUDT15", "rs746071566", "GAGTCG(3)", "GAGTCG(4)")
        .variation("NUDT15", "rs116855232", "C", "T")
    ;
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("NUDT15");
    testWrapper.testPrintCpicCalls("NUDT15", "*1/*2");
    testWrapper.testLookup("NUDT15", "*1", "*2");

    DrugReport azaReport = testWrapper.getContext().getDrugReport(DataSource.CPIC, "azathioprine");
    assertNotNull(azaReport);
    GuidelineReport azaCpicGuideline = azaReport.getGuidelines().iterator().next();
    List<Genotype> genotypes = Genotype.makeGenotypes(azaCpicGuideline.getRelatedGeneReports());
    assertEquals(1, genotypes.size());

    testWrapper.testMatchedAnnotations("azathioprine", 2);
  }

  @Test
  void testNudt15S3(TestInfo testInfo) throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("NUDT15", "rs116855232", "C", "T")
    ;
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("NUDT15");
    testWrapper.testPrintCpicCalls("NUDT15", "*1/*3");
    testWrapper.testLookup("NUDT15", "*1", "*3");

    testWrapper.testMatchedAnnotations("azathioprine", 2);
    testWrapper.testMatchedAnnotations("mercaptopurine", 2);
    testWrapper.testMatchedAnnotations("thioguanine", 2);
  }


  /* MT-RNR1 */
  @Test
  void testMtrnr1(TestInfo testInfo) throws Exception {
    Path outsideCallPath = TestUtils.createTestFile(testInfo, ".tsv");
    try (PrintWriter writer = new PrintWriter(Files.newBufferedWriter(outsideCallPath))) {
      writer.println("MT-RNR1\t1555A>G");
    }
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("CYP2C19")
        .reference("CYP2C9")
    ;
    testWrapper.execute(outsideCallPath);

    testWrapper.testCalledByMatcher("CYP2C19");
    testWrapper.testCalledByMatcher("CYP2C9");
    testWrapper.testReportable("MT-RNR1");
    testWrapper.testMatchedAnnotations("amikacin", 1);
  }


  /* IFNL3/4 */
  @Test
  void testIfnl3(TestInfo testInfo) throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("IFNL3")
    ;
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("IFNL3");
    testWrapper.testReportable("IFNL3");
    testWrapper.testPrintCpicCalls("IFNL3", "rs12979860 reference (C)/rs12979860 reference (C)");
    testWrapper.testMatchedAnnotations("peginterferon alfa-2a", 0);
    testWrapper.testMatchedAnnotations("peginterferon alfa-2b", 0);
  }


  /**
   * This tests the case when an outside call file contains an entry with both a diplotype and a phenotype which is an
   * error.
   */
  @Test
  void testBadOutsideData(TestInfo testInfo) throws Exception {
    Path outsideCallPath = TestUtils.createTestFile(testInfo, ".tsv");
    try (PrintWriter writer = new PrintWriter(Files.newBufferedWriter(outsideCallPath))) {
      writer.println("CYP2D6\t*1/*2\tfoo");
      writer.println("CYP2D6\t*3/*4");
    }

    assertThrows(BadOutsideCallException.class, () -> {
      PharmCATTestWrapper testWrapper = new PharmCATTestWrapper(testInfo, false);
      testWrapper.getVcfBuilder()
          .reference("CYP2C19");
      testWrapper.execute(outsideCallPath);
    });
  }

  @Test
  void testCyp2d6AlleleWithNoFunction(TestInfo testInfo) throws Exception {
    Path outsideCallPath = TestUtils.createTestFile(testInfo,".tsv");
    try (PrintWriter writer = new PrintWriter(Files.newBufferedWriter(outsideCallPath))) {
      writer.println("CYP2D6\t*1/*XXX");
    }

    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("CYP2C19");
    testWrapper.execute(outsideCallPath);

    GeneReport geneReport = testWrapper.getContext().getGeneReport(DataSource.CPIC, "CYP2D6");
    assertNotNull(geneReport);
    assertEquals(1, geneReport.getRecommendationDiplotypes().size());
    Diplotype diplotype = geneReport.getRecommendationDiplotypes().get(0);
    assertEquals("One normal function allele and one unassigned function allele", diplotype.printFunctionPhrase());
  }


  // has multimatch
  @Test
  void testCyp2d6DoubleCall(TestInfo testInfo) throws Exception {
    Path outsideCallPath = TestUtils.createTestFile(testInfo,".tsv");
    try (PrintWriter writer = new PrintWriter(Files.newBufferedWriter(outsideCallPath))) {
      writer.println("CYP2D6\t*1/*1");
      writer.println("CYP2D6\t*1/*2");
    }

    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("CYP2C19");
    testWrapper.execute(outsideCallPath);

    GeneReport geneReport = testWrapper.getContext().getGeneReport(DataSource.CPIC, "CYP2D6");
    assertNotNull(geneReport);
    assertEquals(2, geneReport.getRecommendationDiplotypes().size());
    Diplotype diplotype = geneReport.getRecommendationDiplotypes().get(0);
    assertEquals("Two normal function alleles", diplotype.printFunctionPhrase());
  }

  @Test
  void testCyp3a4(TestInfo testInfo) throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .variation("CYP3A4", "rs72552799", "T", "T")
        .variation("CYP3A4", "rs2242480", "T", "T")
    ;
    testWrapper.execute(null);
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
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper(testInfo, true, true, false);
    testWrapper.getVcfBuilder()
        .reference("CYP2C19")
        .variation("CYP2C19", "rs367543002", "C", "T")
        .variation("CYP2C19", "rs3758581", "G", "G")
        .missing("CYP2C19", "rs367543003");
    testWrapper.execute(null);
    testWrapper.testCalledByMatcher("CYP2C19");
    testWrapper.testReportable("CYP2C19");
  }


  /**
   * Added to have an example of running in CYP2D6-matching mode and make sure messages are applied
   */
  @Test
  void testCallCyp2d6(TestInfo testInfo) throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper(testInfo, false, true, true);
    testWrapper.getVcfBuilder()
        .reference("CYP2D6")
        .reference("CYP2C19")
        .variation("CYP2C19", "rs3758581", "G", "G");
    Path vcfFile = testWrapper.execute(null);
    testWrapper.testCalledByMatcher("CYP2C19", "CYP2D6");
    testWrapper.testReportable("CYP2C19", "CYP2D6");

    Path reporterOutput = vcfFile.getParent().resolve(PharmCAT.getBaseFilename(vcfFile) + ".report.html");
    Document document = Jsoup.parse(reporterOutput.toFile());
    assertNotNull(document.getElementById("CYP2D6"));
    Elements cyp2d6Section = document.select(".gene.CYP2D6");
    assertEquals(1, cyp2d6Section.size());
    assertEquals(0, cyp2d6Section.get(0).getElementsByClass(MessageHelper.MSG_OUTSIDE_CALL).size());
    assertEquals(1, cyp2d6Section.get(0).getElementsByClass(MessageHelper.MSG_CYP2D6_MODE).size());
  }

  /**
   * Can we use activity scores in outside call files? It should be specified in the column for "phenotype"
   */
  @Test
  void testOutsideActivityScore(TestInfo testInfo) throws Exception {
    Path outDir = TestUtils.getTestOutputDir(testInfo, false);
    Path outsideCallPath = outDir.resolve(TestUtils.getTestName(testInfo) + ".tsv");
    try (PrintWriter writer = new PrintWriter(Files.newBufferedWriter(outsideCallPath))) {
      writer.println("CYP2D6\t\t\t1.25");
    }

    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("CYP2C19");
    testWrapper.execute(outsideCallPath);

    testWrapper.testCalledByMatcher("CYP2C19");
    testWrapper.testPrintCalls(DataSource.CPIC, "CYP2C19", "*38/*38");

    testWrapper.testNotCalledByMatcher("CYP2D6");
    testWrapper.testPrintCalls(DataSource.CPIC, "CYP2D6", "Normal Metabolizer (1.25)");

    testWrapper.testMessageCountForGene(DataSource.CPIC, "CYP2C19", 1);
    testWrapper.testGeneHasMessage(DataSource.CPIC, "CYP2C19", "reference-allele");
    testWrapper.testGeneHasMessage(DataSource.CPIC, "CYP2D6", MessageHelper.MSG_OUTSIDE_CALL);
  }

  /**
   * Can we use phenotype in CYP2D6 outside call files?
   */
  @Test
  void testCyp2d6OutsidePhenotype(TestInfo testInfo) throws Exception {
    Path outDir = TestUtils.getTestOutputDir(testInfo, false);
    Path outsideCallPath = outDir.resolve(TestUtils.getTestName(testInfo) + ".tsv");
    try (PrintWriter writer = new PrintWriter(Files.newBufferedWriter(outsideCallPath))) {
      writer.println("CYP2D6\t\tIntermediate Metabolizer");
    }

    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("CYP2C19");
    Path vcfFile = testWrapper.execute(outsideCallPath);

    testWrapper.testCalledByMatcher("CYP2C19");
    testWrapper.testPrintCalls(DataSource.CPIC, "CYP2C19", "*38/*38");

    testWrapper.testNotCalledByMatcher("CYP2D6");
    testWrapper.testPrintCalls(DataSource.CPIC, "CYP2D6", "Intermediate Metabolizer");

    testWrapper.testMessageCountForGene(DataSource.CPIC, "CYP2C19", 1);
    testWrapper.testGeneHasMessage(DataSource.CPIC, "CYP2C19", "reference-allele");
    testWrapper.testGeneHasMessage(DataSource.CPIC, "CYP2D6", MessageHelper.MSG_OUTSIDE_CALL);

    Path reporterOutput = vcfFile.getParent().resolve(PharmCAT.getBaseFilename(vcfFile) + ".report.html");
    Document document = Jsoup.parse(reporterOutput.toFile());
    assertNotNull(document.getElementById("CYP2D6"));
    Elements cyp2d6Section = document.select(".gene.CYP2D6");
    assertEquals(1, cyp2d6Section.size());
    assertEquals(1, cyp2d6Section.get(0).getElementsByClass(MessageHelper.MSG_OUTSIDE_CALL).size());
    assertEquals(0, cyp2d6Section.get(0).getElementsByClass(MessageHelper.MSG_CYP2D6_MODE).size());
  }


  @Test
  void testWarfarinMissingRs12777823(TestInfo testInfo) throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper(testInfo, false);
    testWrapper.getVcfBuilder()
        .reference("CYP2C9")
        .reference("CYP4F2")
        .reference("VKORC1")
        .missingExtraPosition("CYP2C9", "rs12777823")
    ;
    Path vcfFile = testWrapper.execute(null);

    testWrapper.testCalledByMatcher("CYP2C9");
    testWrapper.testReportable("CYP2C9");

    Path reporterOutput = vcfFile.getParent().resolve(PharmCAT.getBaseFilename(vcfFile) + ".report.html");
    Document document = Jsoup.parse(reporterOutput.toFile());
    assertNotNull(document.getElementById("CYP2C9"));
    Element warfarin = document.getElementById("warfarin-cpic-1-1");
    assertNotNull(warfarin);
    Element missing = warfarin.getElementById("warfarin-cpic-1-1-rs12777823");
    assertNotNull(missing);
    assertTrue(missing.text().contains(Haplotype.UNKNOWN));
  }



  private static String printDiagnostic(GeneReport geneReport) {
    return String.format(
        sf_diplotypesTemplate,
        geneReport.getSourceDiplotypes().toString(),
        geneReport.getRecommendationDiplotypes().toString(),
        geneReport.printDisplayCalls()
    );
  }

  
  private static class PharmCATTestWrapper {
    private final Path m_outputPath;
    private final TestVcfBuilder m_vcfBuilder;
    private final boolean m_findCombinations;
    private final boolean m_callCyp2d6;
    private final boolean m_topCandidatesOnly;
    private ReportContext m_reportContext;

    PharmCATTestWrapper(TestInfo testInfo, boolean allMatches) throws IOException {
      this(testInfo, false, allMatches, false);
    }

    PharmCATTestWrapper(TestInfo testInfo, boolean findCombinations, boolean allMatches, boolean callCyp2d6)
        throws IOException {
      Preconditions.checkNotNull(testInfo);

      m_outputPath = TestUtils.getTestOutputDir(testInfo, false);
      if (!Files.isDirectory(m_outputPath)) {
        Files.createDirectories(m_outputPath);
      }
      m_vcfBuilder = new TestVcfBuilder(testInfo).saveFile();
      m_findCombinations = findCombinations;
      m_callCyp2d6 = callCyp2d6;
      m_topCandidatesOnly = !allMatches;
    }

    ReportContext getContext() {
      return m_reportContext;
    }

    TestVcfBuilder getVcfBuilder() {
      return m_vcfBuilder;
    }

    Path execute(Path outsideCallPath) throws Exception {
      Path vcfFile = m_vcfBuilder.generate();
      PharmCAT pcat = new PharmCAT(new Env(),
          true, vcfFile, m_topCandidatesOnly, m_callCyp2d6, m_findCombinations, true,
          true, null, outsideCallPath,
          true, null, null, m_sources, m_compact, true,
          m_outputPath, null, m_compact, PharmCAT.Mode.TEST
      );
      pcat.execute();
      m_reportContext = pcat.getReportContext();
      return vcfFile;
    }


    private void testMatcher(String gene, String... diplotypes) {
      GeneReport geneReport = getContext().getGeneReport(DataSource.CPIC, gene);
      List<String> dips = geneReport.getSourceDiplotypes().stream()
          .map(Diplotype::printBare)
          .sorted()
          .toList();
      try {
        assertThat(dips, contains(diplotypes));
      } catch (AssertionError ex) {
        System.out.println(printDiagnostic(geneReport));
        throw ex;
      }
    }

    /**
     * Test the "print" calls for a gene that will display in the final report or in the phenotyper. This will check
     * that the call count matches and then check each individual call is present (can be 1 or more).
     *
     * @param gene the gene to get diplotypes for
     * @param calls the expected display of the calls, 1 or more
     */
    private void testPrintCpicCalls(String gene, String... calls) {
      testPrintCalls(DataSource.CPIC, gene, calls);
    }

    private void testPrintCalls(DataSource source, String gene, String... calls) {
      GeneReport geneReport = getContext().getGeneReport(source, gene);
      SortedSet<String> dips = new TreeSet<>(geneReport.printDisplayCalls());
      Arrays.sort(calls);
      try {
        assertThat(dips, contains(calls));
      } catch (AssertionError ex) {
        System.out.println(printDiagnostic(geneReport));
        throw ex;
      }
    }

    private void testInferredCalls(String gene, String... calls) {
      GeneReport geneReport = getContext().getGeneReport(DataSource.CPIC, gene);
      SortedSet<String> dips = new TreeSet<>(geneReport.printDisplayInferredCalls());
      Arrays.sort(calls);
      try {
        assertThat(dips, contains(calls));
      } catch (AssertionError ex) {
        System.out.println(printDiagnostic(geneReport));
        throw ex;
      }
    }

    /**
     * Test the diplotype that will be used for looking up the recommendation. This will mostly match what's printed in
     * displays but will differ for particular genes.
     *
     * @param gene the gene to get diplotypes for
     * @param haplotypes the expected haplotypes names used for calling, specifying one will assume homozygous,
     * otherwise specify two haplotype names
     */
    private void testLookup(String gene, String... haplotypes) {
      Preconditions.checkArgument(haplotypes.length >= 1 && haplotypes.length <= 2,
          "Can only test on 1 or 2 haplotypes");

      Map<String, Integer> lookup = new HashMap<>();
      if (haplotypes.length == 2) {
        if (haplotypes[0].equals(haplotypes[1])) {
          lookup.put(haplotypes[0], 2);
        } else {
          lookup.put(haplotypes[0], 1);
          lookup.put(haplotypes[1], 1);
        }
      } else {
        lookup.put(haplotypes[0], 1);
      }


      GeneReport geneReport = getContext().getGeneReport(DataSource.CPIC, gene);
      assertTrue(geneReport.isReportable(), "Not reportable: " + geneReport.getRecommendationDiplotypes());

      assertTrue(geneReport.getRecommendationDiplotypes().stream()
          .anyMatch(d -> d.makeLookupMap().equals(lookup)),
          "Lookup key " + lookup + " not found in lookup " +
              geneReport.getRecommendationDiplotypes().stream().map(Diplotype::makeLookupMap).toList());
    }

    private void testLookupByActivity(String gene, String activityScore) {
      GeneReport geneReport = getContext().getGeneReport(DataSource.CPIC, gene);
      assertTrue(geneReport.isReportable());
      String foundScores = geneReport.getRecommendationDiplotypes().stream()
          .map(Diplotype::getActivityScore)
          .collect(Collectors.joining("; "));
      assertTrue(geneReport.getRecommendationDiplotypes().stream()
          .allMatch(d -> d.printLookupKeys().equals(activityScore)), "Activity score mismatch, expected " + activityScore + " but got " + foundScores);
    }

    /**
     * Check to see if all the given genes have been called by the matcher
     */
    private void testCalledByMatcher(String... genes) {
      assertTrue(genes != null && genes.length > 0);
      Arrays.stream(genes)
          .forEach(g -> {
            assertTrue(getContext().getGeneReports(g).stream().anyMatch(GeneReport::isCalled),
                g + " is not called");
          });
    }

    private void testReportable(String... genes) {
      assertTrue(genes != null && genes.length > 0);
      Arrays.stream(genes)
          .forEach(g -> {
            assertTrue(getContext().getGeneReports(g).stream().anyMatch(GeneReport::isReportable),
                g + " is not reportable");
          });
    }

    /**
     * Check to see if none of the given genes have been called by the matcher
     */
    private void testNotCalledByMatcher(String... genes) {
      Preconditions.checkArgument(genes != null && genes.length > 0);
      Arrays.stream(genes)
          .forEach(g -> {
            assertTrue(getContext().getGeneReports(g).stream().noneMatch(GeneReport::isCalled),
                g + " is called");
          });
    }

    /**
     * Check to see if there is a matching recommendation for the given drug name.
     * @param drugName a drug name that has recommendations
     * @param expectedCount the number of matching recommendations you expect
     */
    private void testMatchedAnnotations(String drugName, int expectedCount) {
      List<DrugReport> drugReports = getContext().getDrugReports(drugName);
      int numMatched = drugReports.stream()
          .mapToInt(DrugReport::getMatchedAnnotationCount)
          .sum();
      assertEquals(expectedCount, numMatched,
          drugName + " has " + numMatched + " matching recommendation(s) instead of " + expectedCount);
    }

    private void testMatchedAnnotations(String drugName, DataSource source, int expectedCount) {
      DrugReport drugReport = getContext().getDrugReport(source, drugName);
      assertNotNull(drugReport);
      assertEquals(expectedCount, drugReport.getMatchedAnnotationCount(),
          drugName + " has " + drugReport.getMatchedAnnotationCount() + " matching " + source +
              " recommendation(s) instead of " + expectedCount);
    }

    private void testAnyMatchFromSource(String drugName, DataSource source) {
      DrugReport drugReport = getContext().getDrugReport(source, drugName);
      assertNotNull(drugReport);
      assertTrue(drugReport.getGuidelines().stream().anyMatch((g) -> g.getSource() == source && g.isMatched()),
          drugName + " does not have matching recommendation from " + source);
    }

    private void testNoMatchFromSource(String drugName, DataSource source) {
      DrugReport drugReport = getContext().getDrugReport(source, drugName);
      if (drugReport != null) {
        assertTrue(drugReport.getGuidelines().stream().noneMatch(r -> r.getSource() == source && r.isMatched()),
            drugName + " has a matching recommendation from " + source + " and expected none");
      }
    }

    private void testMessageCountForDrug(DataSource source, String drugName, int messageCount) {
      DrugReport drugReport = getContext().getDrugReport(source, drugName);
      assertNotNull(drugReport);
      assertEquals(messageCount, drugReport.getMessages().size(),
          drugName + " expected " + messageCount + " messages and got " + drugReport.getMessages());
    }

    private void testMessageCountForGene(DataSource source, String geneName, int messageCount) {
      GeneReport report = getContext().getGeneReport(source, geneName);
      assertNotNull(report);
      assertEquals(messageCount, report.getMessages().size(),
          geneName + " expected " + messageCount + " messages and got " + report.getMessages());
    }

    private void testGeneHasMessage(DataSource source, String geneName, String... msgNames) {
      GeneReport report = getContext().getGeneReport(source, geneName);
      assertNotNull(report);
      assertTrue(report.getMessages().size() >= msgNames.length);
      List<String> names = report.getMessages().stream()
          .map(MessageAnnotation::getName)
          .toList();
      for (String msgName : msgNames) {
        assertTrue(names.contains(msgName),
            geneName + " is missing expected message with name \"" + msgName + "\" (got " + names + ")");
      }
    }
  }
}

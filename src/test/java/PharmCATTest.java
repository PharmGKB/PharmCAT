import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import org.apache.commons.io.FileUtils;
import org.jsoup.Jsoup;
import org.jsoup.nodes.Document;
import org.jsoup.nodes.Element;
import org.jsoup.select.Elements;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.TestInfo;
import org.pharmgkb.common.util.PathUtils;
import org.pharmgkb.pharmcat.ParseException;
import org.pharmgkb.pharmcat.PharmCAT;
import org.pharmgkb.pharmcat.ReportableException;
import org.pharmgkb.pharmcat.TestUtils;
import org.pharmgkb.pharmcat.TestVcfBuilder;
import org.pharmgkb.pharmcat.haplotype.ResultSerializer;
import org.pharmgkb.pharmcat.haplotype.model.GeneCall;
import org.pharmgkb.pharmcat.haplotype.model.Result;
import org.pharmgkb.pharmcat.phenotype.Phenotyper;
import org.pharmgkb.pharmcat.reporter.ReportContext;
import org.pharmgkb.pharmcat.reporter.model.DataSource;
import org.pharmgkb.pharmcat.reporter.model.MessageAnnotation;
import org.pharmgkb.pharmcat.reporter.model.VariantReport;
import org.pharmgkb.pharmcat.reporter.model.result.Diplotype;
import org.pharmgkb.pharmcat.reporter.model.result.DrugReport;
import org.pharmgkb.pharmcat.reporter.model.result.GeneReport;

import static com.github.stefanbirkner.systemlambda.SystemLambda.tapSystemOut;
import static org.junit.jupiter.api.Assertions.*;


/**
 * Test the data generated from a full run of the PharmCAT matcher and reporter.
 *
 * @author Ryan Whaley
 */
class PharmCATTest {
  private static final Path sf_outputDir = TestUtils.TEST_OUTPUT_DIR.resolve(PharmCATTest.class.getSimpleName());
  private static final boolean SAVE_TEST_OUTPUT = "true".equals(System.getenv("PHARMCAT_SAVE_TEST_OUTPUT"));
  private static final String sf_outsideCalls = """
      ##Test Outside Call Data
      #Gene\tDiplotype\tdiplotype activity\tdiplotype calling notes\tjaccard\tpart\tpValue\tROI notes\tspecial case\tnomenclature version
      CYP2D6\tCYP2D6*1/CYP2D6*4\t\t\t0.6\t0.75\tp: 0.0\t\t\tv1.9-2017_02_09
      """;
  private static final String sf_otherOutsideCalls = "CYP2D6\t*3/*4\nG6PD\tB (wildtype)/B (wildtype)\n";
  private static final String sf_diplotypesTemplate = "\nmatcher: %s\nreporter: %s\nprint (displayCalls): %s";
  private static Path s_outsideCallFilePath;
  private static Path s_otherOutsideCallFilePath;


  @BeforeAll
  static void prepare() throws IOException {

    s_outsideCallFilePath = TestUtils.createTempFile("outsideCall", ".tsv");
    try (FileWriter fw = new FileWriter(s_outsideCallFilePath.toFile())) {
      fw.write(sf_outsideCalls);
    }

    s_otherOutsideCallFilePath = TestUtils.createTempFile("otherOutsideCall", ".tsv");
    try (FileWriter fw = new FileWriter(s_otherOutsideCallFilePath.toFile())) {
      fw.write(sf_otherOutsideCalls);
    }
  }

  @AfterAll
  static void teardown() {
    if (!SAVE_TEST_OUTPUT) {
      try {
        FileUtils.deleteDirectory(sf_outputDir.toFile());
      } catch (IOException ex) {
        // ignore
      }
    }
  }


  @Test
  void cliBasic() throws Exception {
    Path vcfFile = PathUtils.getPathToResource("reference.vcf");
    Path refMatcherOutput = vcfFile.getParent().resolve("reference.match.json");
    Path refPhenotyperOutput = vcfFile.getParent().resolve("reference.phenotype.json");
    Path refReporterOutput = vcfFile.getParent().resolve("reference.report.html");

    // require VCF
    try {
      String systemOut = tapSystemOut(() -> {
        PharmCAT.main(null);
      });
      //System.out.println(systemOut);
      assertTrue(systemOut.contains("No input"));
      assertTrue(systemOut.contains("-vcf"));
    } finally {
      Files.deleteIfExists(refMatcherOutput);
      Files.deleteIfExists(refPhenotyperOutput);
      Files.deleteIfExists(refReporterOutput);
    }

    // standard full run
    try {
      String systemOut = tapSystemOut(() -> {
        PharmCAT.main(new String[] {"-vcf", vcfFile.toString()});
      });
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

      List<GeneReport> reports = Phenotyper.readGeneReports(refPhenotyperOutput);
      Optional<GeneReport> grOpt = reports.stream()
          .filter(gr -> gr.getGene().equals("CYP2D6"))
          .findFirst();
      assertTrue(grOpt.isPresent());
      assertFalse(grOpt.get().isCalled());

      Document document = Jsoup.parse(refReporterOutput.toFile());
      assertNotNull(document.getElementById("genotypes"));
      assertNotNull(document.getElementById("CYP2D6"));
      assertNotNull(document.getElementById("PA166104937"));

    } finally {
      Files.deleteIfExists(refMatcherOutput);
      Files.deleteIfExists(refPhenotyperOutput);
      Files.deleteIfExists(refReporterOutput);
    }
  }


  @Test
  void cliCallCyp2d6(TestInfo testInfo) throws Exception {
    Path vcfFile = PathUtils.getPathToResource("reference.vcf");

    Path outputDir = TestUtils.getTestOutputDir(testInfo, true);
    Path refMatcherOutput = outputDir.resolve("reference.match.json");
    Path refPhenotyperOutput = outputDir.resolve("reference.phenotype.json");
    Path refReporterOutput = outputDir.resolve("reference.report.html");

    try {
      String systemOut = tapSystemOut(() -> {
        PharmCAT.main(new String[]{
            "-vcf", vcfFile.toString(),
            "-o", outputDir.toString(),
            "-research", "cyp2d6" });
      });
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

      List<GeneReport> reports = Phenotyper.readGeneReports(refPhenotyperOutput);
      Optional<GeneReport> grOpt = reports.stream()
          .filter(gr -> gr.getGene().equals("CYP2D6"))
          .findFirst();
      assertTrue(grOpt.isPresent());
      assertTrue(grOpt.get().isCalled());
      assertFalse(grOpt.get().isOutsideCall());

      Document document = Jsoup.parse(refReporterOutput.toFile());
      assertNotNull(document.getElementById("genotypes"));
      Element geneTitle = document.getElementById("CYP2D6");
      assertNotNull(geneTitle);
      Elements diplotypes = geneTitle.parent().getElementsByTag("li");
      for (String diplotype : diplotypes.eachText()) {
        assertEquals("*1/*1", diplotype);
      }
      assertNotNull(document.getElementById("PA166104937"));

    } finally {
      Files.deleteIfExists(refMatcherOutput);
      Files.deleteIfExists(refPhenotyperOutput);
      Files.deleteIfExists(refReporterOutput);
    }
  }

  @Test
  void cliCallCyp2d6WithOverlappingOutsideCall(TestInfo testInfo) throws Exception {
    Path vcfFile = PathUtils.getPathToResource("reference.vcf");
    Path outsideCallFile = PathUtils.getPathToResource("PharmCATTest-cyp2d6.tsv");

    Path outputDir = TestUtils.getTestOutputDir(testInfo, true);
    Path refMatcherOutput = outputDir.resolve("reference.match.json");
    Path refPhenotyperOutput = outputDir.resolve("reference.phenotype.json");
    Path refReporterOutput = outputDir.resolve("reference.report.html");

    try {
      String systemOut = tapSystemOut(() -> {
        PharmCAT.main(new String[]{
            "-vcf", vcfFile.toString(),
            "-po", outsideCallFile.toString(),
            "-o", outputDir.toString(),
            "-research", "cyp2d6" });
      });
      System.out.println(systemOut);
      assertTrue(systemOut.contains("Cannot specify outside call for CYP2D6, it is already called in sample data"));
      assertTrue(Files.exists(refMatcherOutput));
      assertFalse(Files.exists(refPhenotyperOutput));
      assertFalse(Files.exists(refReporterOutput));

    } finally {
      Files.deleteIfExists(refMatcherOutput);
      Files.deleteIfExists(refPhenotyperOutput);
      Files.deleteIfExists(refReporterOutput);
    }
  }

  @Test
  void cliOutsideCallsOnly(TestInfo testInfo) throws Exception {
    Path outsideCallFile = PathUtils.getPathToResource("PharmCATTest-cyp2d6.tsv");

    Path outputDir = TestUtils.getTestOutputDir(testInfo, true);
    Path outsideMatcherOutput = outputDir.resolve("PharmCATTest-cyp2d6.match.json");
    Path outsidePhenotyperOutput = outputDir.resolve("PharmCATTest-cyp2d6.phenotype.json");
    Path outsideReporterOutput = outputDir.resolve("PharmCATTest-cyp2d6.report.html");

    try {
      String systemOut = tapSystemOut(() -> {
        PharmCAT.main(new String[] {
            "-o", outputDir.toString(),
            "-phenotyper",
            "-po", outsideCallFile.toString()
        });
      });
      System.out.println(systemOut);
      assertTrue(systemOut.contains("Done."));
      assertFalse(Files.exists(outsideMatcherOutput));
      assertTrue(Files.exists(outsidePhenotyperOutput));
      assertFalse(Files.exists(outsideReporterOutput));

      List<GeneReport> reports = Phenotyper.readGeneReports(outsidePhenotyperOutput);
      Optional<GeneReport> grOpt = reports.stream()
          .filter(gr -> gr.getGene().equals("CYP2D6"))
          .findFirst();
      assertTrue(grOpt.isPresent());
      assertFalse(grOpt.get().isCalled());
      assertTrue(grOpt.get().isOutsideCall());

    } finally {
      Files.deleteIfExists(outsideMatcherOutput);
      Files.deleteIfExists(outsidePhenotyperOutput);
      Files.deleteIfExists(outsideReporterOutput);
    }
  }

  /**
   * This test is dependent on CYP2C19 definition.
   * It may fail if CYP2C19 definition changes too much.
   */
  @Test
  void cliMatchAllFlag(TestInfo testInfo) throws Exception {
    Path vcfFile = PathUtils.getPathToResource("PharmCATTest-cyp2c19MissingPositions.vcf");

    Path outputDir = TestUtils.getTestOutputDir(testInfo, true);
    String baseFilename = TestUtils.getTestName(testInfo);
    Path matcherOutput = outputDir.resolve(baseFilename + ".match.json");
    Path phenotyperOutput = outputDir.resolve(baseFilename + ".phenotype.json");
    Path reporterOutput = outputDir.resolve(baseFilename + ".report.html");

    // matcher only, expecting 1 CYP2C19 matches
    try {
      String systemOut = tapSystemOut(() -> {
        PharmCAT.main(new String[] {
            "-vcf", vcfFile.toString(),
            "-matcher",
            "-o", outputDir.toString(),
            "-bf", baseFilename,
        });
      });
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
      Files.deleteIfExists(matcherOutput);
      Files.deleteIfExists(phenotyperOutput);
      Files.deleteIfExists(reporterOutput);
    }

    // matcher only, expecting many CYP2C19 matches
    try {
      String systemOut = tapSystemOut(() -> {
        PharmCAT.main(new String[] {
            "-vcf", vcfFile.toString(),
            "-matcher",
            "-ma",
            "-o", outputDir.toString(),
            "-bf", baseFilename,
        });
      });
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
      Files.deleteIfExists(matcherOutput);
      Files.deleteIfExists(phenotyperOutput);
      Files.deleteIfExists(reporterOutput);
    }
  }


  /**
   * This test ensures that the output of the phenotyper is consistent regardless of whether you run the matcher first
   * in one process then the phenotyper in a second process versus if you run the matcher and phenotyper together in
   * one process.
   */
  @Test
  void cliTestConsistentOutput() throws Exception {
    Path vcfFile = PathUtils.getPathToResource("reference.vcf");
    Path singlesMatcherOutput = vcfFile.getParent().resolve("singles.match.json");
    Path singlesPhenotyperOutput = vcfFile.getParent().resolve("singles.phenotype.json");
    Path doublePhenotyperOutput = vcfFile.getParent().resolve("double.phenotype.json");

    try {
      String systemOut = tapSystemOut(() -> {
        PharmCAT.main(new String[] {"-matcher", "-bf", "singles", "-vcf", vcfFile.toString()});
      });
      assertTrue(systemOut.contains("Done."));
      assertTrue(Files.exists(singlesMatcherOutput));

      String singlesPhenoOut = tapSystemOut(() -> {
        PharmCAT.main(new String[]{"-phenotyper", "-pi", singlesMatcherOutput.toString(), "-bf", "singles"});
      });
      assertTrue(singlesPhenoOut.contains("Done."));
      assertTrue(Files.exists(singlesPhenotyperOutput));

      String doubleOut = tapSystemOut(() -> {
        PharmCAT.main(new String[]{
            "-matcher",
            "-phenotyper",
            "-vcf", vcfFile.toString(),
            "-bf", "double"});
      });
      assertTrue(doubleOut.contains("Done."));
      assertTrue(Files.exists(doublePhenotyperOutput));

      assertEquals(
          Files.readString(singlesPhenotyperOutput),
          Files.readString(doublePhenotyperOutput));
    } finally {
      Files.deleteIfExists(singlesMatcherOutput);
      Files.deleteIfExists(singlesPhenotyperOutput);
      Files.deleteIfExists(doublePhenotyperOutput);
    }
  }


  /**
   * NOTE: if these assertions fail then new data may have been added from the DataManager because of an update to the
   * CPIC database. If that's true, then update these numbers to the current count. If the count changes with no known
   * change to the CPIC database then something may be wrong in code.
   */
  @Test
  void testCounts() throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper("top", false);
    testWrapper.getVcfBuilder()
            .reference("CYP2C9");
    testWrapper.execute(null);
    assertEquals(23, testWrapper.getContext().getGeneReports().size());
    assertEquals(93, testWrapper.getContext().getDrugReports().size());
  }

  @Test
  void testAll() throws Exception {
    Path outsideCallPath = TestUtils.createTempFile("hlab", ".tsv");
    try (FileWriter fw = new FileWriter(outsideCallPath.toFile())) {
      fw.write(
          "CYP2D6\t*3/*4\n" +
          "G6PD\tB (wildtype)/B (wildtype)\n" +
          "HLA-A\t\t*31:01 positive\n" +
          "HLA-B\t*15:02/*57:01"
      );
    }
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper("all", false);
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
   * able to come up with a matched group
   */
  @Test
  void testCyp2c19() throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper("cyp2c19.singleGeneMatch", false);
    testWrapper.getVcfBuilder()
        .variation("CYP2C19", "rs3758581", "G", "G");
    testWrapper.execute(s_otherOutsideCallFilePath);

    testWrapper.testCalledByMatcher("CYP2C19");
    testWrapper.testPrintCalls( "CYP2C19", "*1/*1");

    testWrapper.testMatchedGroups("amitriptyline", 2);
    testWrapper.testAnyMatchFromSource("amitriptyline", DataSource.CPIC);
    testWrapper.testAnyMatchFromSource("amitriptyline", DataSource.DPWG);
    testWrapper.testMatchedGroups("citalopram", 2);
    testWrapper.testAnyMatchFromSource("citalopram", DataSource.CPIC);
    testWrapper.testAnyMatchFromSource("citalopram", DataSource.DPWG);
    testWrapper.testMatchedGroups("ivacaftor", 0);
  }

  /**
   * Tests how PharmCAT handles that state when sample VCF data exists for a gene and an outside call also exists for
   * that gene. Currently, this should execute successfully by ignoring outside call data and using the sample data
   */
  @Test
  void testCallerCollision() throws Exception {
    Path outsideCallPath = Files.createTempFile("cyp2c19_collision", ".tsv");
    try (FileWriter fw = new FileWriter(outsideCallPath.toFile())) {
      fw.write("CYP2C19\t*2/*2\n");
    }

    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper("caller.collision", false);
    testWrapper.getVcfBuilder()
        .reference("CYP2C19");
    testWrapper.execute(outsideCallPath);

    testWrapper.testCalledByMatcher("CYP2C19");
    // this is the diplotype indicated in the VCF, not the one in the outside call
    testWrapper.testPrintCalls( "CYP2C19", "*38/*38");

    GeneReport geneReport = testWrapper.getContext().getGeneReport("CYP2C19");
    assertEquals(1, geneReport.getMessages().size());
    assertTrue(geneReport.getMessages().stream().allMatch(m -> m.getExceptionType().equals(MessageAnnotation.TYPE_NOTE)));
  }

  /**
   * This test case demos that an "ambiguity" {@link MessageAnnotation} which specifies a variant and a diplotype call
   * for a given drug report will be matched and added to the {@link DrugReport}
   */
  @Test
  void testCyp2c19_s1s2rs58973490het() throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper("cyp2c19.s1s2rs58973490het", false);
    testWrapper.getVcfBuilder()
        .variation("CYP2C19", "rs12769205", "A", "G")
        .variation("CYP2C19", "rs58973490", "G", "A")
        .variation("CYP2C19", "rs4244285", "G", "A")
        .variation("CYP2C19", "rs3758581", "G", "G");
    testWrapper.execute(s_otherOutsideCallFilePath);

    testWrapper.testCalledByMatcher("CYP2C19");
    testWrapper.testPrintCalls( "CYP2C19", "*1/*2");
    testWrapper.testNotCalledByMatcher("CYP2D6");
    testWrapper.testPrintCalls( "CYP2D6", "*3/*4");

    testWrapper.testMatchedGroups("amitriptyline", 2);
    testWrapper.testAnyMatchFromSource("amitriptyline", DataSource.CPIC);
    testWrapper.testAnyMatchFromSource("amitriptyline", DataSource.DPWG);
    testWrapper.testMatchedGroups("citalopram", 2);
    testWrapper.testAnyMatchFromSource("citalopram", DataSource.CPIC);
    testWrapper.testAnyMatchFromSource("citalopram", DataSource.DPWG);
    testWrapper.testMatchedGroups("clomipramine", 3);
    testWrapper.testAnyMatchFromSource("clomipramine", DataSource.CPIC);
    testWrapper.testAnyMatchFromSource("clomipramine", DataSource.DPWG);
    testWrapper.testMatchedGroups("ivacaftor", 0);

    GeneReport cyp2c19report = testWrapper.getContext().getGeneReport("CYP2C19");

    VariantReport vr = cyp2c19report.findVariantReport("rs58973490")
        .orElseThrow(() -> new RuntimeException("Variant missing from test data"));
    assertTrue(vr.isHetCall());

    // ambiguity message will not apply in this case because all variants are available for CYP2C19, but one message
    // should appear for the *1 call
    assertEquals(1, cyp2c19report.getMessages().stream()
        .filter(m -> m.getExceptionType().equals(MessageAnnotation.TYPE_AMBIGUITY) && m.getMatches().getVariant().equals("rs58973490"))
        .count());

    testWrapper.testMessageCountForDrug("amitriptyline", 2);
  }

  /**
   * This test case demos that an "ambiguity" {@link MessageAnnotation} which specifies a variant and a diplotype call
   * for a given drug report will not be matched when the variant in the message is homozygous
   */
  @Test
  void testCyp2c19_s1s2() throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper("cyp2c19.s1s2", false);
    testWrapper.getVcfBuilder()
        .variation("CYP2C19", "rs12769205", "A", "G")
        .variation("CYP2C19", "rs4244285", "G", "A")
        .variation("CYP2C19", "rs3758581", "G", "G");
    testWrapper.execute(s_otherOutsideCallFilePath);

    testWrapper.testCalledByMatcher("CYP2C19");
    testWrapper.testPrintCalls( "CYP2C19", "*1/*2");

    testWrapper.testMatchedGroups("amitriptyline", 2);
    testWrapper.testAnyMatchFromSource("amitriptyline", DataSource.CPIC);
    testWrapper.testAnyMatchFromSource("amitriptyline", DataSource.DPWG);
    testWrapper.testMatchedGroups("citalopram", 2);
    testWrapper.testAnyMatchFromSource("citalopram", DataSource.CPIC);
    testWrapper.testAnyMatchFromSource("citalopram", DataSource.DPWG);
    testWrapper.testMatchedGroups("clomipramine", 3);
    testWrapper.testAnyMatchFromSource("clomipramine", DataSource.CPIC);
    testWrapper.testAnyMatchFromSource("clomipramine", DataSource.DPWG);
    testWrapper.testMatchedGroups("ivacaftor", 0);

    GeneReport cyp2c19report = testWrapper.getContext().getGeneReport("CYP2C19");

    // make sure the variant in question is not a het call
    VariantReport vr = cyp2c19report.findVariantReport("rs58973490")
        .orElseThrow(() -> new RuntimeException("Variant missing from test data"));
    assertFalse(vr.isHetCall());

    // the variant is hom so ambiguity message should not apply and, thus, no matching messages
    assertEquals(0, cyp2c19report.getMessages().stream()
        .filter(m -> m.getExceptionType().equals(MessageAnnotation.TYPE_AMBIGUITY) && m.getMatches().getVariant().equals("rs58973490"))
        .count());

    DrugReport amiReport = testWrapper.getContext().getDrugReport("amitriptyline");
    // should only get the *1 message, the variant is hom so ambiguity message should not match
    assertEquals(1, amiReport.getMessages().size());
  }

  @Test
  void testClomipramineCall() throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper("cyp2c19.clomipramine", false);
    testWrapper.getVcfBuilder()
        .variation("CYP2C19", "rs12769205", "G", "G")
        .variation("CYP2C19", "rs4244285", "A", "A")
        .variation("CYP2C19", "rs3758581", "G", "G");
    testWrapper.execute(s_otherOutsideCallFilePath);

    testWrapper.testCalledByMatcher("CYP2C19");
    testWrapper.testPrintCalls( "CYP2C19", "*2/*2");

    testWrapper.testMatchedGroups("amitriptyline", 2);
    testWrapper.testAnyMatchFromSource("amitriptyline", DataSource.CPIC);
    testWrapper.testAnyMatchFromSource("amitriptyline", DataSource.DPWG);
    testWrapper.testMatchedGroups("clomipramine", 3);
    testWrapper.testAnyMatchFromSource("clomipramine", DataSource.CPIC);
    testWrapper.testAnyMatchFromSource("clomipramine", DataSource.DPWG);
    testWrapper.testMatchedGroups("desipramine", 1);
    testWrapper.testAnyMatchFromSource("desipramine", DataSource.CPIC);
    testWrapper.testMatchedGroups("doxepin", 2);
    testWrapper.testAnyMatchFromSource("doxepin", DataSource.CPIC);
    testWrapper.testAnyMatchFromSource("doxepin", DataSource.DPWG);
    testWrapper.testMatchedGroups("imipramine", 3);
    testWrapper.testAnyMatchFromSource("imipramine", DataSource.CPIC);
    testWrapper.testAnyMatchFromSource("imipramine", DataSource.DPWG);
    testWrapper.testMatchedGroups("nortriptyline", 2);
    testWrapper.testAnyMatchFromSource("nortriptyline", DataSource.CPIC);
    testWrapper.testAnyMatchFromSource("nortriptyline", DataSource.DPWG);
    testWrapper.testMatchedGroups("trimipramine", 1);

    testWrapper.testMatchedGroups("clopidogrel", 4);
    testWrapper.testAnyMatchFromSource("clopidogrel", DataSource.CPIC);
    testWrapper.testAnyMatchFromSource("clopidogrel", DataSource.DPWG);

    testWrapper.testMatchedGroups("lansoprazole", 2);
    testWrapper.testAnyMatchFromSource("lansoprazole", DataSource.CPIC);
    testWrapper.testAnyMatchFromSource("lansoprazole", DataSource.DPWG);

    // voriconazole has 2 populations with recommendations so should have 2 matching groups from CPIC, and 1 from DPWG
    testWrapper.testMatchedGroups("voriconazole", 3);
    testWrapper.testAnyMatchFromSource("voriconazole", DataSource.CPIC);
    testWrapper.testAnyMatchFromSource("voriconazole", DataSource.DPWG);
  }

  @Test
  void testCyp2c19noCall() throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper("cyp2c19.noCall", false);
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
  void testCyp2c19s4bs17rs28399504missing() throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper("cyp2c19.s4bs17rs28399504missing", false);
    testWrapper.getVcfBuilder()
        .variation("CYP2C19", "rs12248560", "T", "T")
        .missing("CYP2C19", "rs28399504")
        .variation("CYP2C19", "rs3758581", "G", "G");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("CYP2C19");
    testWrapper.testPrintCalls("CYP2C19", "*4/*4", "*4/*17", "*17/*17");

    testWrapper.testMatchedGroups("citalopram", 6);
    testWrapper.testMatchedGroups("citalopram", DataSource.CPIC, 3);
    testWrapper.testMatchedGroups("citalopram", DataSource.DPWG, 3);
    testWrapper.testAnyMatchFromSource("citalopram", DataSource.CPIC);
    testWrapper.testAnyMatchFromSource("citalopram", DataSource.DPWG);
  }

  @Test
  void testCyp2c19s1s4het() throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper("cyp2c19.s1s4het", false);
    testWrapper.getVcfBuilder()
        .variation("CYP2C19", "rs12248560", "T", "T")
        .variation("CYP2C19", "rs28399504", "A", "G")
        .variation("CYP2C19", "rs3758581", "G", "G");
    testWrapper.execute(s_outsideCallFilePath);

    testWrapper.testCalledByMatcher("CYP2C19");
    testWrapper.testReportable("CYP2D6");

    testWrapper.testPrintCalls("CYP2D6", "*1/*4");
    testWrapper.testPrintCalls("CYP2C19", "*4/*17");

    assertTrue(testWrapper.getContext().getGeneReport("CYP2D6").isOutsideCall());
  }

  @Test
  void testCyp2c19s1s4missingS1() throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper("cyp2c19.s1s4missingS1", false);
    testWrapper.getVcfBuilder()
        .variation("CYP2C19", "rs12248560", "C", "T")
        .variation("CYP2C19", "rs28399504", "A", "G")
        .missing("CYP2C19", "rs3758581");
    testWrapper.execute(s_outsideCallFilePath);

    testWrapper.testCalledByMatcher("CYP2C19");
    testWrapper.testReportable("CYP2D6");

    testWrapper.testPrintCalls("CYP2D6", "*1/*4");
    testWrapper.testPrintCalls("CYP2C19",  "*1/*4", "*4/*38");

    assertTrue(testWrapper.getContext().getGeneReport("CYP2D6").isOutsideCall());
    GeneReport cyp2c19report = testWrapper.getContext().getGeneReport("CYP2C19");
    assertEquals("Yes", cyp2c19report.isMissingVariants());

    assertFalse(cyp2c19report.isPhased());
    assertTrue(cyp2c19report.findVariantReport("rs12248560").map(VariantReport::isHetCall).orElse(false));
    assertTrue(cyp2c19report.findVariantReport("rs3758581").map(VariantReport::isMissing).orElse(false));

    assertTrue(cyp2c19report.hasHaplotype("*38"));

    // message is for *1/*4 being ambiuguous with unphased data
    assertEquals(2, cyp2c19report.getMessages().stream()
        .filter(m -> m.getExceptionType().equals(MessageAnnotation.TYPE_AMBIGUITY))
        .count());

    DrugReport amiReport = testWrapper.getContext().getDrugReport("amitriptyline");
    assertEquals(3, amiReport.getMessages().size());
  }

  @Test
  void testCyp2c19s4s17() throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper("cyp2c19.s1s4missingS1", false);
    testWrapper.getVcfBuilder()
        .variation("CYP2C19", "rs12248560", "C", "T")
        .variation("CYP2C19", "rs28399504", "A", "G")
        .variation("CYP2C19", "rs3758581", "G", "G");
    testWrapper.execute(s_outsideCallFilePath);

    testWrapper.testCalledByMatcher("CYP2C19");
    testWrapper.testReportable("CYP2D6");

    testWrapper.testPrintCalls("CYP2D6", "*1/*4");
    testWrapper.testPrintCalls("CYP2C19", "*1/*4");

    assertTrue(testWrapper.getContext().getGeneReport("CYP2D6").isOutsideCall());
  }

  @Test
  void testCftrRefRef() throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper("cftr.ref_ref", false);
    testWrapper.getVcfBuilder()
        .reference("CFTR");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("CFTR");
    testWrapper.testPrintCalls("CFTR", "No CPIC variants found");
    testWrapper.testLookup("CFTR", "ivacaftor non-responsive CFTR sequence");

    testWrapper.testMatchedGroups("ivacaftor", 1);
  }

  @Test
  void testCftrD1270NHet() throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper("cftr.ref_D1270N", false);
    testWrapper.getVcfBuilder()
        .variation("CFTR", "rs11971167", "G", "A");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("CFTR");
    testWrapper.testPrintCalls("CFTR", "D1270N (heterozygous)");
    testWrapper.testLookup("CFTR", "ivacaftor non-responsive CFTR sequence", "D1270N");

    testWrapper.testMatchedGroups("ivacaftor", 1);
  }

  @Test
  void testCftrD1270NG551D() throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper("cftr.ref_D1270NG551D", false);
    testWrapper.getVcfBuilder()
        .variation("CFTR", "rs11971167", "G", "A")
        .variation("CFTR", "rs75527207", "G", "A");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("CFTR");
    testWrapper.testPrintCalls("CFTR", "D1270N/G551D");
    testWrapper.testLookup("CFTR", "G551D", "D1270N");

    testWrapper.testMatchedGroups("ivacaftor", 1);
  }

  @Test
  void testRosuvastatin() throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper("rosuvastatin", false);
    testWrapper.getVcfBuilder()
        .variation("ABCG2", "rs2231142", "G", "T")
        .variation("SLCO1B1", "rs56101265", "T", "C");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("ABCG2", "SLCO1B1");
    testWrapper.testPrintCalls("SLCO1B1", "*1/*2");

    testWrapper.testMatchedGroups("rosuvastatin", 1);
  }

  @Test
  void testAmitryptylineCallWoCyp2c19() throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper("AmitryptylineCallWoCyp2c19", false);
    testWrapper.getVcfBuilder()
        .reference("DPYD");
    testWrapper.execute(s_outsideCallFilePath);

    testWrapper.testReportable("CYP2D6");
    testWrapper.testPrintCalls("CYP2D6", "*1/*4");
    testWrapper.testLookup("CYP2D6", "*1", "*4");

    testWrapper.testMatchedGroups("amitriptyline", 2);
    testWrapper.testAnyMatchFromSource("amitriptyline", DataSource.CPIC);
    testWrapper.testAnyMatchFromSource("amitriptyline", DataSource.DPWG);
  }

  @Test
  void testSlco1b1HomWild() throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper("slco1b1.s1s1", false);
    testWrapper.getVcfBuilder()
        .reference("SLCO1B1");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("SLCO1B1");
    testWrapper.testPrintCalls("SLCO1B1", "*1/*1");
    testWrapper.testLookup("SLCO1B1", "*1");

    GeneReport slco1b1Report = testWrapper.getContext().getGeneReport("SLCO1B1");
    assertTrue(slco1b1Report.getHighlightedVariants().contains("rs4149056T/rs4149056T"));
  }

  @Test
  void testSlco1b1HomVar() throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper("slco1b1.s5s15", false);
    testWrapper.getVcfBuilder()
        .variation("SLCO1B1", "rs2306283", "A", "G")
        .variation("SLCO1B1", "rs4149056", "C", "C");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("SLCO1B1");
    testWrapper.testPrintCalls("SLCO1B1", "*5/*15");
    testWrapper.testLookup("SLCO1B1", "*5", "*15");
  }

  @Test
  void testSlco1b1Test5() throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper("slco1b1.s1s44", false);
    testWrapper.getVcfBuilder()
        .variation("SLCO1B1", "rs2306283", "A", "G")
        .variation("SLCO1B1", "rs11045852", "A", "G")
        .variation("SLCO1B1", "rs74064213", "A", "G");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("SLCO1B1");
    testWrapper.testPrintCalls("SLCO1B1", "*1/*44");
    testWrapper.testLookup("SLCO1B1", "*1", "*44");
  }

  @Test
  void testSlco1b1Test3() throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper("slco1b1.s1s15", false);
    testWrapper.getVcfBuilder()
        .variation("SLCO1B1", "rs2306283", "A", "G")
        .variation("SLCO1B1", "rs4149056", "T", "C");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("SLCO1B1");
    testWrapper.testPrintCalls("SLCO1B1", "*1/*15");
    testWrapper.testLookup("SLCO1B1", "*1", "*15");
  }

  @Test
  void testSlco1b1Test4() throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper("slco1b1.s5s45", false);
    testWrapper.getVcfBuilder()
        .variation("SLCO1B1", "rs4149056", "T", "C")
        .variation("SLCO1B1", "rs71581941", "C", "T");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("SLCO1B1");
    testWrapper.testPrintCalls("SLCO1B1", "*5/*45");
    testWrapper.testLookup("SLCO1B1", "*5", "*45");

    testWrapper.testMatchedGroups("simvastatin", 1);
  }

  @Test
  void testDpydS1S2B() throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper("dpyd.s1s2b", false);
    testWrapper.getVcfBuilder()
        .phased()
        .variation("DPYD", "rs3918290", "C", "T")
        .variation("DPYD", "rs1801159", "C", "T");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("DPYD");

    GeneReport dpyd = testWrapper.getContext().getGeneReport("DPYD");
    assertEquals(2, dpyd.getMatcherAlleles().size());
    testWrapper.testPrintCalls("DPYD", "c.1627A>G (*5)/c.1905+1G>A (*2A)");
    testWrapper.testLookup("DPYD", "c.1627A>G (*5)", "c.1905+1G>A (*2A)");

    testWrapper.testMatchedGroups("fluorouracil", 1);
    testWrapper.testMatchedGroups("capecitabine", 1);
  }

  @Test
  void testDpydUnphased() throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper("dpyd.unphased", false);
    testWrapper.getVcfBuilder()
        .variation("DPYD", "rs3918290", "C", "T")
        .variation("DPYD", "rs1801159", "C", "T");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("DPYD");

    GeneReport dpyd = testWrapper.getContext().getGeneReport("DPYD");
    assertEquals(3, dpyd.getMatcherAlleles().size());
    assertEquals(1, dpyd.getReporterDiplotypes().size());
    testWrapper.testPrintCalls("DPYD", "c.1627A>G (*5)", "c.1905+1G>A (*2A)", "Reference");
    testWrapper.testLookupByActivity("DPYD", "1");

    testWrapper.testMatchedGroups("fluorouracil", 1);
    testWrapper.testMatchedGroups("capecitabine", 1);
  }

  @Test
  void testDpydUnphasedMultiple() throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper("dpyd.unphased.multiple", false);
    testWrapper.getVcfBuilder()
        .variation("DPYD", "rs183385770", "C", "T")  // 0 activity value
        .variation("DPYD", "rs186169810", "A", "C") // 0.5 activity value
        .variation("DPYD", "rs112766203", "G", "A") // 0.5 activity value
        .variation("DPYD", "rs144395748", "G", "C"); // 1 activity value
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("DPYD");

    GeneReport dpyd = testWrapper.getContext().getGeneReport("DPYD");
    assertEquals(5, dpyd.getMatcherAlleles().size());
    assertEquals(1, dpyd.getReporterDiplotypes().size());
    testWrapper.testPrintCalls("DPYD", "Reference", "c.1024G>A", "c.1314T>G", "c.1358C>G", "c.2279C>T");
    testWrapper.testLookupByActivity("DPYD", "0.5");

    testWrapper.testMatchedGroups("fluorouracil", 1);
    testWrapper.testMatchedGroups("capecitabine", 1);
  }

  @Test
  void testDpydC2846het() throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper("dpyd.c2846het", false);
    testWrapper.getVcfBuilder()
        .phased()
        .variation("DPYD", "rs67376798", "T", "A");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("DPYD");
    testWrapper.testPrintCalls("DPYD", "Reference/c.2846A>T");
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
  void testDpydPhasedMultiTrans() throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper("dpyd.phased.multi.trans", false);
    testWrapper.getVcfBuilder()
        .phased()
        .variation("DPYD", "rs67376798", "A", "T") // Strand 1 decreased - c.2846A>T
        .variation("DPYD", "rs72547601", "C", "T") // Strand 1 no function - c.2933A>G
        .variation("DPYD", "rs60139309", "T", "C") // Strand 2 normal function - c.2582A>G
        .variation("DPYD", "rs139834141", "C", "T") // Strand 2 normal function - c.498G>A
    ;
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("DPYD");
    testWrapper.testPrintCalls("DPYD", "c.498G>A + c.2582A>G/c.2846A>T + c.2933A>G");
    testWrapper.testLookup("DPYD", "c.2933A>G", "c.498G>A");

    testWrapper.testAnyMatchFromSource("fluorouracil", DataSource.CPIC);

    testWrapper.testAnyMatchFromSource("capecitabine", DataSource.CPIC);
  }

  /**
   * This test is the same as the previous test but DPYD is unphased instead of phased. This means the individual found
   * alleles should be reported and then the two least-function alleles should be used for recommendation lookup.
   */
  @Test
  void testDpydUnphasedMulti() throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper("dpyd.unphased.multi", false);
    testWrapper.getVcfBuilder()
        .variation("DPYD", "rs67376798", "A", "T") // decreased - c.2846A>T
        .variation("DPYD", "rs72547601", "C", "T") // no function - c.2933A>G
        .variation("DPYD", "rs60139309", "T", "C") // normal function - c.2582A>G
        .variation("DPYD", "rs139834141", "C", "T") // normal function - c.498G>A
    ;
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("DPYD");
    testWrapper.testPrintCalls("DPYD", "c.498G>A", "c.2582A>G", "c.2846A>T", "c.2933A>G", "Reference");
    testWrapper.testLookup("DPYD", "c.2933A>G", "c.2846A>T");

    testWrapper.testAnyMatchFromSource("fluorouracil", DataSource.CPIC);

    testWrapper.testAnyMatchFromSource("capecitabine", DataSource.CPIC);
  }

  @Test
  void testDpydS12het() throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper("dpyd.c1156het", false);
    testWrapper.getVcfBuilder()
        .phased()
        .variation("DPYD", "rs78060119", "C", "A");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("DPYD");
    testWrapper.testPrintCalls("DPYD", "Reference/c.1156G>T (*12)");
    testWrapper.testLookup("DPYD", "Reference", "c.1156G>T (*12)");

    testWrapper.testMatchedGroups("fluorouracil", 1);
    testWrapper.testMatchedGroups("capecitabine", 1);
  }

  /**
   * This tests a special case of SLCO1B1. The gene in this scenario is "uncalled" by the matcher due the sample VCF
   * data. However, SLCO1B1 has an override that will display the rs4149056 diplotype regardless of call status. That
   * same override will assign alleles to use for recommendation lookup
   */
  @Test
  void testSlco1b1TestMulti() throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper("slco1b1.multi", false);
    testWrapper.getVcfBuilder()
        .variation("SLCO1B1", "rs2306283", "G", "G")
        .variation("SLCO1B1", "rs4149056", "T", "C")
        .variation("SLCO1B1", "rs11045853", "A", "A")
        .variation("SLCO1B1", "rs72559748", "G", "G");
    testWrapper.execute(null);

    testWrapper.testNotCalledByMatcher("SLCO1B1");
    testWrapper.testPrintCalls("SLCO1B1", "rs4149056T/rs4149056C");
    testWrapper.testLookup("SLCO1B1", "*1", "*5");

    testWrapper.testMatchedGroups("simvastatin", 1);
  }

  @Test
  void testUgt1a1Phased() throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper("ugt1a1.s1s80.phased", false);
    testWrapper.getVcfBuilder()
        .phased()
        .variation("UGT1A1", "rs887829", "C", "T");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("UGT1A1");
    testWrapper.testPrintCalls("UGT1A1", "*1/*80");
    testWrapper.testLookup("UGT1A1", "*1", "*80");
  }

  @Test
  void testUgt1a1Unphased() throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper("ugt1a1.s1s80.unphased", false);
    testWrapper.getVcfBuilder()
        .phased()
        .variation("UGT1A1", "rs887829", "C", "T");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("UGT1A1");
    testWrapper.testPrintCalls("UGT1A1", "*1/*80");
    testWrapper.testLookup("UGT1A1", "*1", "*80");
  }

  @Test
  void testUgt1a1s1s1() throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper("ugt1a1.s1s1", false);
    testWrapper.getVcfBuilder()
        .reference("UGT1A1");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("UGT1A1");
    testWrapper.testPrintCalls("UGT1A1", "*1/*1");
    testWrapper.testLookup("UGT1A1", "*1", "*1");
  }

  @Test
  void testUgt1a1S1S80S28() throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper("ugt1a1.s1s80s28unphased", false);
    testWrapper.getVcfBuilder()
        .variation("UGT1A1", "rs887829", "C", "T")
        .variation("UGT1A1", "rs3064744", "TA(7)", "TA(8)");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("UGT1A1");
    testWrapper.testPrintCalls("UGT1A1", "*1/*80+*28");
    testWrapper.testLookup("UGT1A1", "*1", "*80+*28");
  }

  @Test
  void testUgt1a1S28S37() throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper("ugt1a1.s28s37", false);
    testWrapper.getVcfBuilder()
        .variation("UGT1A1", "rs3064744", "TA(8)", "TA(9)");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("UGT1A1");
    testWrapper.testPrintCalls("UGT1A1", "*28/*37");
    testWrapper.testLookup("UGT1A1", "*37", "*28");
  }

  @Test
  void testUgt1a1s28s80phased() throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper("ugt1a1.s1s80s28phased", false);
    testWrapper.getVcfBuilder()
        .phased()
        .variation("UGT1A1", "rs887829", "C", "T")
        .variation("UGT1A1", "rs3064744", "TA(7)", "TA(8)");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("UGT1A1");
    testWrapper.testPrintCalls("UGT1A1", "*1/*80+*28");
    testWrapper.testLookup("UGT1A1", "*1", "*80+*28");

    // the guideline should have a matching message for the *1 call but no ambiguity call
    testWrapper.testMatchedGroups("atazanavir", 1);
    testWrapper.testMessageCountForDrug("atazanavir", 1);
  }

  @Test
  void testUgt1a1s28s80s6s60phased() throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper("ugt1a1.s28s80s6s60phased", false);
    testWrapper.getVcfBuilder()
        .phased()
        .variation("UGT1A1", "rs887829", "T", "C")
        .variation("UGT1A1", "rs3064744", "TA(8)", "TA(7)")
        .variation("UGT1A1", "rs4148323", "G", "A");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("UGT1A1");
    testWrapper.testPrintCalls("UGT1A1", "*6/*80+*28");
    testWrapper.testLookup("UGT1A1", "*6", "*80+*28");
  }

  @Test
  void testUgt1a1s28s80s6s60unphased() throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper("ugt1a1.s28s80s6s60unphased", false);
    testWrapper.getVcfBuilder()
        .variation("UGT1A1", "rs887829", "T", "C")
        .variation("UGT1A1", "rs3064744", "TA(8)", "TA(7)")
        .variation("UGT1A1", "rs4148323", "G", "A");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("UGT1A1");
    testWrapper.testPrintCalls("UGT1A1", "*6/*80+*28");
    testWrapper.testLookup("UGT1A1", "*6", "*80+*28");
  }

  @Test
  void testUgt1a1s6s6() throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper("ugt1a1.s6s6", false);
    testWrapper.getVcfBuilder()
        .variation("UGT1A1", "rs4148323", "A", "A");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("UGT1A1");
    testWrapper.testPrintCalls("UGT1A1", "*6/*6");
    testWrapper.testLookup("UGT1A1", "*6", "*6");
  }

  @Test
  void testUgt1a1s6s60s80s28MissingPhased() throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper("ugt1a1.s6s60s80s28MissingPhased", false);
    testWrapper.getVcfBuilder()
        .phased()
        .missing("UGT1A1", "rs3064744")
        .variation("UGT1A1", "rs887829", "C", "T")
        .variation("UGT1A1", "rs4148323", "A", "G");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("UGT1A1");
    testWrapper.testPrintCalls("UGT1A1", "*6/*80", "*6/*80+*28", "*6/*80+*37");
    testWrapper.testLookup("UGT1A1", "*6", "*80");
    testWrapper.testLookup("UGT1A1", "*6", "*80+*28");
    testWrapper.testLookup("UGT1A1", "*6", "*80+*37");
  }

  @Test
  void testUgt1a1s6s60s80s28MissingUnphased() throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper("ugt1a1.s6s60s80s28MissingUnphased", false);
    testWrapper.getVcfBuilder()
        .missing("UGT1A1", "rs3064744")
        .variation("UGT1A1", "rs887829", "C", "T")
        .variation("UGT1A1", "rs4148323", "A", "G");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("UGT1A1");
    testWrapper.testPrintCalls("UGT1A1", "*6/*80", "*6/*80+*28", "*6/*80+*37");
    testWrapper.testLookup("UGT1A1", "*6", "*80");
    testWrapper.testLookup("UGT1A1", "*6", "*80+*28");
    testWrapper.testLookup("UGT1A1", "*6", "*80+*37");
  }

  @Test
  void testUgt1a1s80s28missing() throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper("ugt1a1.s80s28missing", false);
    testWrapper.getVcfBuilder()
        .missing("UGT1A1", "rs3064744")
        .variation("UGT1A1", "rs887829", "C", "T");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("UGT1A1");
    testWrapper.testPrintCalls("UGT1A1", "*1/*80", "*1/*80+*28", "*1/*80+*37");
    testWrapper.testLookup("UGT1A1", "*1", "*80");
    testWrapper.testLookup("UGT1A1", "*1", "*80+*28");
    testWrapper.testLookup("UGT1A1", "*1", "*80+*37");
  }

  @Test
  void testUgt1a1na12717() throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper("ugt1a1.na12717", false);
    testWrapper.getVcfBuilder()
        .variation("UGT1A1", "rs887829", "T", "T")
        .variation("UGT1A1", "rs3064744", "TA(8)", "TA(7)");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("UGT1A1");
    testWrapper.testPrintCalls("UGT1A1", "*80/*80+*28");
    testWrapper.testLookup("UGT1A1", "*80", "*80+*28");
  }

  @Test
  void testUgt1a1s28homMissing() throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper("ugt1a1.s28s28unphaseds60s80miss", false);
    testWrapper.getVcfBuilder()
        .missing("UGT1A1", "rs887829")
        .variation("UGT1A1", "rs3064744", "TA(8)", "TA(8)");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("UGT1A1");
    testWrapper.testPrintCalls("UGT1A1", "*28/*28", "*28/*80+*28", "*80+*28/*80+*28");
    testWrapper.testLookup("UGT1A1", "*28", "*28");
    testWrapper.testLookup("UGT1A1", "*28", "*80+*28");
    testWrapper.testLookup("UGT1A1", "*80+*28", "*80+*28");
  }

  @Test
  void testUgt1a1s28s60Hom() throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper("ugt1a1.s1s28hom", false);
    testWrapper.getVcfBuilder()
        .variation("UGT1A1", "rs3064744", "TA(8)", "TA(7)");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("UGT1A1");
    testWrapper.testPrintCalls("UGT1A1", "*1/*28");
    testWrapper.testLookup("UGT1A1", "*1", "*28");
  }

  @Test
  void testUgt1a1s27s28unphaseds80s60missing() throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper("ugt1a1.s27s28unphaseds80s60missing", false);
    testWrapper.getVcfBuilder()
        .missing("UGT1A1", "rs887829")
        .variation("UGT1A1", "rs3064744", "TA(8)", "TA(7)")
        .variation("UGT1A1", "rs35350960", "C", "A");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("UGT1A1");
    testWrapper.testPrintCalls("UGT1A1", "*27/*28", "*27/*80+*28");
    testWrapper.testLookup("UGT1A1", "*27", "*28");
    testWrapper.testLookup("UGT1A1", "*27", "*80+*28");
  }

  @Test
  void testUgt1a1HG00436() throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper("ugt1a1.HG00436", false);
    testWrapper.getVcfBuilder()
        .phased()
        .variation("UGT1A1", "rs887829", "T", "C")
        .variation("UGT1A1", "rs3064744", "TA(8)", "TA(7)")
        .variation("UGT1A1", "rs35350960", "A", "C");
    testWrapper.execute(null);

    testWrapper.testNotCalledByMatcher("UGT1A1");
  }

  @Test
  void testUgt1a1s1s80s27s60s28missingphased() throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper("ugt1a1.s1s80s27s60s28missingphased", false);
    testWrapper.getVcfBuilder()
        .phased()
        .missing("UGT1A1", "rs3064744")
        .variation("UGT1A1", "rs887829", "T", "C")
        .variation("UGT1A1", "rs35350960", "A", "C");
    testWrapper.execute(null);

    testWrapper.testNotCalledByMatcher("UGT1A1");
    GeneReport geneReport = testWrapper.getContext().getGeneReport("UGT1A1");
    assertTrue(geneReport.isPhased());
  }

  @Test
  void testUgt1a1s1s60s80s6phased() throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper("ugt1a1.s1s60s80s6phased", false);
    testWrapper.getVcfBuilder()
        .phased()
        .variation("UGT1A1", "rs887829", "T", "C")
        .variation("UGT1A1", "rs35350960", "A", "C");
    testWrapper.execute(null);

    testWrapper.testNotCalledByMatcher("UGT1A1");
    GeneReport geneReport = testWrapper.getContext().getGeneReport("UGT1A1");
    assertTrue(geneReport.isPhased());
  }

  @Test
  void testUgt1a1s1s60s80s28s6phased() throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper("ugt1a1.s1s60s80s28s6phased", false);
    testWrapper.getVcfBuilder()
        .phased()
        .variation("UGT1A1", "rs887829", "T", "C")
        .variation("UGT1A1", "rs3064744", "TA(8)", "TA(7)")
        .variation("UGT1A1", "rs35350960", "A", "C");
    testWrapper.execute(null);

    testWrapper.testNotCalledByMatcher("UGT1A1");
    GeneReport geneReport = testWrapper.getContext().getGeneReport("UGT1A1");
    assertTrue(geneReport.isPhased());
  }

  @Test
  void testUgt1a1s1s37s80s60phased() throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper("ugt1a1.s1s37s80s60phased", false);
    testWrapper.getVcfBuilder()
        .phased()
        .variation("UGT1A1", "rs887829", "T", "C")
        .variation("UGT1A1", "rs3064744", "TA(9)", "TA(7)");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("UGT1A1");
    testWrapper.testReportable("UGT1A1");
    testWrapper.testPrintCalls("UGT1A1", "*1/*80+*37");
    testWrapper.testLookup("UGT1A1", "*1", "*80+*37");
    GeneReport geneReport = testWrapper.getContext().getGeneReport("UGT1A1");
    assertTrue(geneReport.isPhased());
  }

  @Test
  void testCyp3a5Missing3Message() throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper("cyp3a5.s3missing", false);
    testWrapper.getVcfBuilder()
        .missing("CYP3A5", "rs776746");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("CYP3A5");
    testWrapper.testReportable("CYP3A5");
    testWrapper.testPrintCalls("CYP3A5", "*1/*1");
    testWrapper.testLookup("CYP3A5", "*1", "*1");

    GeneReport gene = testWrapper.getContext().getGeneReport("CYP3A5");
    // rs776746 should be missing from this report
    assertNotNull(gene.getVariantReports());
    assertTrue(gene.getVariantReports().stream().anyMatch(v -> v.isMissing() && v.getDbSnpId().equals("rs776746")));

    // the guideline should have a matching message
    assertTrue(testWrapper.getContext().getDrugReports().stream()
        .filter(r -> r.getRelatedDrugs().contains("tacrolimus"))
        .allMatch(r -> r.getMessages().size() > 0));

    assertFalse(gene.isPhased());
  }

  @Test
  void testCyp3a5v1() throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper("cyp3a5.s1s3rs776746rs55965422het", false);
    testWrapper.getVcfBuilder()
        .variation("CYP3A5", "rs776746", "T", "C");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("CYP3A5");
    testWrapper.testReportable("CYP3A5");
    testWrapper.testPrintCalls("CYP3A5", "*1/*3");
    testWrapper.testLookup("CYP3A5", "*1", "*3");
  }

  @Test
  void testCyp3a5v2() throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper("cyp3a5.v2", false);
    testWrapper.getVcfBuilder()
        .variation("CYP3A5", "rs28383479", "C", "T")
        .variation("CYP3A5", "rs776746", "C", "T")
    ;
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("CYP3A5");
    testWrapper.testPrintCalls("CYP3A5", "*3/*9");
    testWrapper.testLookup("CYP3A5", "*3", "*9");
  }

  @Test
  void testCyp3a5v3() throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper("cyp3a5.v3", false);
    testWrapper.getVcfBuilder()
        .variation("CYP3A5", "rs776746", "C", "C")
    ;
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("CYP3A5");
    testWrapper.testPrintCalls("CYP3A5", "*3/*3");
    testWrapper.testLookup("CYP3A5", "*3", "*3");
  }

  @Test
  void testCyp3a5v4() throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper("cyp3a5.v4", false);
    testWrapper.getVcfBuilder()
        .variation("CYP3A5", "rs776746", "T", "C")
    ;
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("CYP3A5");
    testWrapper.testPrintCalls("CYP3A5", "*1/*3");
    testWrapper.testLookup("CYP3A5", "*1", "*3");
  }

  @Test
  void testCyp3a5v5() throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper("cyp3a5.v4", false);
    testWrapper.getVcfBuilder()
        .variation("CYP3A5", "rs28383479", "T", "C")
        .variation("CYP3A5", "rs776746", "T", "C")
    ;
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("CYP3A5");
    testWrapper.testPrintCalls("CYP3A5", "*3/*9");
    testWrapper.testLookup("CYP3A5", "*3", "*9");
  }

  @Test
  void testHlab(TestInfo testInfo) throws Exception {
    Path outsideCallPath = TestUtils.createTempFile(testInfo, ".tsv");
    try (FileWriter fw = new FileWriter(outsideCallPath.toFile())) {
      fw.write("HLA-B\t*15:02/*57:01");
    }

    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper("hlab", false);
    testWrapper.getVcfBuilder()
        .reference("CYP2C9");
    testWrapper.execute(outsideCallPath);

    testWrapper.testCalledByMatcher("CYP2C9");
    testWrapper.testNotCalledByMatcher("HLA-B");
    testWrapper.testReportable("CYP2C9");
    testWrapper.testReportable("HLA-B");
    testWrapper.testMatchedGroups("abacavir", DataSource.CPIC, 1);
    testWrapper.testMatchedGroups("abacavir", DataSource.DPWG, 1);
    testWrapper.testMatchedGroups("allopurinol", DataSource.CPIC, 1);
    testWrapper.testMatchedGroups("allopurinol", DataSource.DPWG, 1);
    testWrapper.testMatchedGroups("phenytoin", 4);
    testWrapper.testAnyMatchFromSource("phenytoin", DataSource.CPIC);
    testWrapper.testAnyMatchFromSource("phenytoin", DataSource.DPWG);
  }

  @Test
  void testHlabPhenotype(TestInfo testInfo) throws Exception {
    Path outsideCallPath = TestUtils.createTempFile(testInfo, ".tsv");
    try (FileWriter fw = new FileWriter(outsideCallPath.toFile())) {
      fw.write("HLA-B\t\t*57:01 positive");
    }

    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper("hlab.phenotype", false);
    testWrapper.getVcfBuilder()
        .reference("CYP2C9");
    testWrapper.execute(outsideCallPath);

    testWrapper.testCalledByMatcher("CYP2C9");
    testWrapper.testNotCalledByMatcher("HLA-B");
    testWrapper.testReportable("CYP2C9");
    testWrapper.testReportable("HLA-B");
    testWrapper.testMatchedGroups("abacavir", 2);
    testWrapper.testMatchedGroups("abacavir", DataSource.CPIC, 1);
    testWrapper.testMatchedGroups("abacavir", DataSource.DPWG, 1);
    // allopurinol relies on a different allele for recs so no matches
    testWrapper.testMatchedGroups("allopurinol", 0);
    // phenytoin also relies on a different allele but there will be a match for DPWG since the recommendations are
    // split between the two genes on that side
    testWrapper.testMatchedGroups("phenytoin", 1);
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
  void testRecommendationExamples() throws Exception {
    Path outsideCallPath = Files.createTempFile("noFunction", ".tsv");
    try (FileWriter fw = new FileWriter(outsideCallPath.toFile())) {
      fw.write("HLA-A\t\t*31:01 positive\n");
      fw.write("HLA-B\t*57:01/*58:01\t\n");
    }

    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper("recommendation.examples", false);
    testWrapper.getVcfBuilder()
        .reference("CYP2C9")
        .variation("CYP2C19", "rs12769205", "G", "G")
        .variation("CYP2C19", "rs4244285", "A", "A")
        .variation("CYP2C19", "rs3758581", "G", "G");
    testWrapper.execute(outsideCallPath);

    testWrapper.testLookup("CYP2C19", "*2", "*2");
    testWrapper.testPrintCalls("CYP2C19", "*2/*2");

    GeneReport cyp2c9 = testWrapper.getContext().getGeneReport("CYP2C9");
    assertEquals(1, cyp2c9.getReporterDiplotypes().size());
    assertTrue(cyp2c9.getReporterDiplotypes().stream().allMatch(d -> d.getActivityScore().equals("2.0")));

    testWrapper.testReportable("CYP2C19", "CYP2C9", "HLA-A", "HLA-B");
    testWrapper.testMatchedGroups("celecoxib", 1);
    testWrapper.testAnyMatchFromSource("celecoxib", DataSource.CPIC);
    testWrapper.testMatchedGroups("citalopram", 2);
    testWrapper.testMatchedGroups("clomipramine", 2);
    testWrapper.testMatchedGroups("clopidogrel", 4);
    testWrapper.testMatchedGroups("clopidogrel", DataSource.CPIC, 3);
    testWrapper.testMatchedGroups("clopidogrel", DataSource.DPWG, 1);
    testWrapper.testNoMatchFromSource("flucloxacillin", DataSource.CPIC);
    testWrapper.testMatchedGroups("flucloxacillin", DataSource.DPWG, 1);
    testWrapper.testNoMatchFromSource("fluvoxamine", DataSource.CPIC);
    testWrapper.testNoMatchFromSource("fluvoxamine", DataSource.DPWG);
    testWrapper.testMatchedGroups("siponimod", 1);
    testWrapper.testAnyMatchFromSource("siponimod", DataSource.DPWG);

    testWrapper.testMatchedGroups("carbamazepine", 5);
  }

  @Test
  void testTpmtStar1s() throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper("tpmt.star1s", false);
    testWrapper.getVcfBuilder()
        .variation("TPMT", "rs1800460", "C", "T")
        .variation("TPMT", "rs1142345", "T", "C");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("TPMT");
    testWrapper.testPrintCalls("TPMT", "*1/*3A");
    testWrapper.testLookup("TPMT", "*1", "*3A");

    GeneReport tpmtReport = testWrapper.getContext().getGeneReport("TPMT");
    assertEquals(43, tpmtReport.getVariantReports().size());
    assertEquals(0, tpmtReport.getHighlightedVariants().size());
  }


  @Test
  void testCyp2c9star61() throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper("cyp2c9.s1s61", false);
    testWrapper.getVcfBuilder()
        .variation("CYP2C9", "rs1799853", "C", "T")
        .variation("CYP2C9", "rs202201137", "A", "G");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("CYP2C9");
    testWrapper.testPrintCalls("CYP2C9", "*1/*61");
    testWrapper.testLookup("CYP2C9", "*1", "*61");
  }

  @Test
  void testCyp2c9star1Hom() throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper("cyp2c9.s1s1", false);
    testWrapper.getVcfBuilder()
        .reference("CYP2C9");
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("CYP2C9");
    testWrapper.testPrintCalls("CYP2C9", "*1/*1");
    testWrapper.testLookup("CYP2C9", "*1", "*1");
    testWrapper.testMatchedGroups("celecoxib", 1);
    testWrapper.testMatchedGroups("ibuprofen", 1);
    testWrapper.testMatchedGroups("lornoxicam", 1);
  }


  /**
   * Test CYP2B6 for a het *34 sample file. When doing the "top match" scenario this will only match to a 1/34 and,
   * thus, only match to a single recommendation. This test will have a different outcome when run in "all matches" mode
   * and should be compared with {@link #testCyp2b6star1star34AllMatch()}.
   */
  @Test
  void testCyp2b6star1star34() throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper("cyp2b6.s1s34", false);
    testWrapper.getVcfBuilder()
        .variation("CYP2B6", "rs34223104", "T", "C")
        .variation("CYP2B6", "rs3211371", "C", "A")
        .variation("CYP2B6", "rs3745274", "G", "T")
        .variation("CYP2B6", "rs2279343", "A", "G")
    ;
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("CYP2B6");
    testWrapper.testPrintCalls("CYP2B6", "*1/*34");
    testWrapper.testLookup("CYP2B6", "*1", "*34");
    testWrapper.testMatchedGroups("efavirenz", 1);
  }

  /**
   * This test is just like {@link #testCyp2b6star1star34()} but run in "all matches" mode. This should result in 2
   * possible different calls coming from the matcher. These two have different phenotypes and, thus, match to different
   * recommendations.
   */
  @Test
  void testCyp2b6star1star34AllMatch() throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper("cyp2b6.s1s34.allMatches", true);
    testWrapper.getVcfBuilder()
        .variation("CYP2B6", "rs34223104", "T", "C")
        .variation("CYP2B6", "rs3211371", "C", "A")
        .variation("CYP2B6", "rs3745274", "G", "T")
        .variation("CYP2B6", "rs2279343", "A", "G")
    ;
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("CYP2B6");
    testWrapper.testPrintCalls("CYP2B6", "*1/*34", "*33/*36");
    testWrapper.testLookup("CYP2B6", "*1", "*34");
    testWrapper.testLookup("CYP2B6", "*33", "*36");
    testWrapper.testMatchedGroups("efavirenz", 2);
  }


  /* NUDT15 */
  @Test
  void testNudt15Ref() throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper("nudt15.s1s1", false);
    testWrapper.getVcfBuilder()
        .reference("NUDT15");
    testWrapper.execute(null);

    testWrapper.testPrintCalls("NUDT15", "*1/*1");
    testWrapper.testLookup("NUDT15", "*1", "*1");

    testWrapper.testMatchedGroups("azathioprine", 2);
    testWrapper.testMatchedGroups("mercaptopurine", 2);
    testWrapper.testAnyMatchFromSource("mercaptopurine", DataSource.CPIC);
    testWrapper.testAnyMatchFromSource("mercaptopurine", DataSource.DPWG);
    testWrapper.testMatchedGroups("thioguanine", 2);
    testWrapper.testAnyMatchFromSource("thioguanine", DataSource.CPIC);
    testWrapper.testAnyMatchFromSource("thioguanine", DataSource.DPWG);
  }

  @Test
  void testNudt15S2() throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper("nudt15.s1s2", false);
    testWrapper.getVcfBuilder()
        .variation("NUDT15", "rs746071566", "GAGTCG(3)", "GAGTCG(4)")
        .variation("NUDT15", "rs116855232", "C", "T")
    ;
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("NUDT15");
    testWrapper.testPrintCalls("NUDT15", "*1/*2");
    testWrapper.testLookup("NUDT15", "*1", "*2");

    testWrapper.testMatchedGroups("azathioprine", 2);
  }

  @Test
  void testNudt15S3() throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper("nudt15.s1s3", false);
    testWrapper.getVcfBuilder()
        .variation("NUDT15", "rs116855232", "C", "T")
    ;
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("NUDT15");
    testWrapper.testPrintCalls("NUDT15", "*1/*3");
    testWrapper.testLookup("NUDT15", "*1", "*3");

    testWrapper.testMatchedGroups("azathioprine", 2);
    testWrapper.testMatchedGroups("mercaptopurine", 2);
    testWrapper.testMatchedGroups("thioguanine", 2);
  }


  /* MT-RNR1 */
  @Test
  void testMtrnr1() throws Exception {
    Path outsideCallPath = Files.createTempFile("mtrnr1", ".tsv");
    try (FileWriter fw = new FileWriter(outsideCallPath.toFile())) {
      fw.write("MT-RNR1\t1555A>G\n");
    }
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper("mtrnr1", false);
    testWrapper.getVcfBuilder()
        .reference("CYP2C19")
        .reference("CYP2C9")
    ;
    testWrapper.execute(outsideCallPath);

    testWrapper.testCalledByMatcher("CYP2C19");
    testWrapper.testCalledByMatcher("CYP2C9");
    testWrapper.testReportable("MT-RNR1");
    testWrapper.testMatchedGroups("amikacin", 1);
  }


  /* IFNL3/4 */
  @Test
  void testIfnl3() throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper("ifnl3", false);
    testWrapper.getVcfBuilder()
        .reference("IFNL3")
    ;
    testWrapper.execute(null);

    testWrapper.testCalledByMatcher("IFNL3");
    testWrapper.testReportable("IFNL3");
    testWrapper.testPrintCalls("IFNL3", "rs12979860 reference (C)/rs12979860 reference (C)");
    testWrapper.testMatchedGroups("peginterferon alfa-2a", 0);
    testWrapper.testMatchedGroups("peginterferon alfa-2b", 0);
  }


  /**
   * This tests the case when an outside call file contains an entry with both a diplotype and a phenotype which is an
   * error.
   */
  @Test
  void testBadOutsideData(TestInfo testInfo) throws Exception {
    Path badOutsideDataPath = TestUtils.createTempFile(testInfo, ".tsv");
    try (FileWriter fw = new FileWriter(badOutsideDataPath.toFile())) {
      fw.write("CYP2D6\t*1/*2\tfoo\nCYP2D6\t*3/*4");
    }

    try {
      PharmCATTestWrapper testWrapper = new PharmCATTestWrapper("bad.data", false);
      testWrapper.getVcfBuilder()
              .reference("CYP2C19");
      testWrapper.execute(badOutsideDataPath);
      fail("Should have failed due to a duplicate gene definition in outside call");
    }
    catch (ParseException ex) {
      // we want this to fail so ignore handling the exception
    }
  }

  @Test
  void testCyp2d6AlleleWithNoFunction(TestInfo testInfo) throws Exception {
    Path outsideCallPath = TestUtils.createTempFile(testInfo,".tsv");
    try (FileWriter fw = new FileWriter(outsideCallPath.toFile())) {
      fw.write("CYP2D6\t*1/*XXX\n");
    }

    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper("bad.data", false);
    testWrapper.getVcfBuilder()
        .reference("CYP2C19");
    testWrapper.execute(outsideCallPath);

    GeneReport geneReport = testWrapper.getContext().getGeneReport("CYP2D6");
    assertNotNull(geneReport);
    assertEquals(1, geneReport.getReporterDiplotypes().size());
    Diplotype diplotype = geneReport.getReporterDiplotypes().get(0);
    assertEquals("One normal function allele and one unassigned function allele", diplotype.printFunctionPhrase());
  }

  @Test
  void testCyp2d6DoubleCall(TestInfo testInfo) throws Exception {
    Path outsideCallPath = TestUtils.createTempFile(testInfo,".tsv");
    try (FileWriter fw = new FileWriter(outsideCallPath.toFile())) {
      fw.write("CYP2D6\t*1/*1\nCYP2D6\t*1/*2\n");
    }

    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper("double.outside.call", false);
    testWrapper.getVcfBuilder()
        .reference("CYP2C19");
    testWrapper.execute(outsideCallPath);

    GeneReport geneReport = testWrapper.getContext().getGeneReport("CYP2D6");
    assertNotNull(geneReport);
    assertEquals(2, geneReport.getReporterDiplotypes().size());
    Diplotype diplotype = geneReport.getReporterDiplotypes().get(0);
    assertEquals("Two normal function alleles", diplotype.printFunctionPhrase());
  }

  /**
   * This tests the case when the incoming sample file has coverage for a gene but then an outside call is also
   * supplied. This should result in an error
   */
  @Test
  void testCallerCollision(TestInfo testInfo) throws Exception {
    Path outsidePath = TestUtils.createTempFile(testInfo, ".tsv");
    try (FileWriter fw = new FileWriter(outsidePath.toFile())) {
      fw.write("CYP2C19\t*2/*2");
    }

    try {
      PharmCATTestWrapper testWrapper = new PharmCATTestWrapper("conflict.two.data.sources", false);
      testWrapper.getVcfBuilder()
          .reference("CYP2C19");
      testWrapper.execute(outsidePath);
      fail("Should have failed due to a duplicate gene definition between matcher and outside caller");
    }
    catch (ReportableException ex) {
      // we want this to fail so ignore handling the exception
    }
  }

  @Test
  void testCyp3a4() throws Exception {
    PharmCATTestWrapper testWrapper = new PharmCATTestWrapper("cyp3a4", false);
    testWrapper.getVcfBuilder()
        .variation("CYP3A4", "rs72552799", "T", "T")
        .variation("CYP3A4", "rs2242480", "T", "T")
    ;
    testWrapper.execute(null);
    testWrapper.testCalledByMatcher("CYP3A4");
    testWrapper.testReportable("CYP3A4");
    testWrapper.testPrintCalls("CYP3A4", "*8/*8");
    testWrapper.testMatchedGroups("quetiapine", 1);
  }


  private static String printDiagnostic(GeneReport geneReport) {
    return String.format(
        sf_diplotypesTemplate,
        geneReport.getMatcherDiplotypes().toString(),
        geneReport.getReporterDiplotypes().toString(),
        geneReport.printDisplayCalls()
    );
  }

  
  private static class PharmCATTestWrapper {
    private final PharmCAT m_pharmcat;
    private final Path m_outputPath;
    private final TestVcfBuilder m_vcfBuilder;

    PharmCATTestWrapper(String testKey, boolean allMatches) throws IOException {
      this(testKey, false, allMatches);
    }

    PharmCATTestWrapper(String testKey, boolean findCombinations, boolean allMatches) throws IOException {

      m_outputPath = sf_outputDir.resolve(testKey);
      if (!Files.isDirectory(m_outputPath)) {
        Files.createDirectories(m_outputPath);
      }
      m_vcfBuilder = new TestVcfBuilder(testKey).saveFile();

      m_pharmcat = new PharmCAT(true)
          .matchCombinations(findCombinations)
          .matchTopCandidateOnly(!allMatches);
    }

    ReportContext getContext() {
      return m_pharmcat.getReportContext();
    }

    TestVcfBuilder getVcfBuilder() {
      return m_vcfBuilder;
    }

    void execute(Path outsideCallPath) throws Exception {
      Path vcfFile = m_vcfBuilder.generate(m_outputPath);
      m_pharmcat.execute(vcfFile, outsideCallPath);
    }

    /**
     * Test the "print" calls for a gene that will display in the final report or in the phenotyper. This will check that
     * the call count matches and then check each individual call is present (can be 1 or more).
     * @param gene the gene to get diplotypes for
     * @param calls the expected display of the calls, 1 or more
     */
    private void testPrintCalls(String gene, String... calls) {
      GeneReport geneReport = getContext().getGeneReport(gene);
      Collection<String> dips = geneReport.printDisplayCalls();
      assertEquals(calls.length, dips.size(), "Expected " + gene + " call count (" + calls.length + ") doesn't match actual call count (" + dips.size() + "): " + String.join(", ", dips));
      Arrays.stream(calls)
          .forEach(c -> assertTrue(dips.contains(c),
              c + " not in " + gene + ":" + dips + printDiagnostic(geneReport)));
    }

    /**
     * Test the diplotype that will be used for looking up the recommendation. This will mostly match what's printed in
     * displays but will differ for particular genes
     * @param gene the gene to get diplotypes for
     * @param haplotypes the expected haplotypes names used for calling, specifying one will assume homozygous, otherwise specify two haplotype names
     */
    private void testLookup(String gene, String... haplotypes) {
      Map<String, Integer> lookup = new HashMap<>();
      if (haplotypes.length == 1) {
        lookup.put(haplotypes[0], 2);
      } else if (haplotypes.length == 2) {
        if (haplotypes[0].equals(haplotypes[1])) {
          lookup.put(haplotypes[0], 2);
        } else {
          lookup.put(haplotypes[0], 1);
          lookup.put(haplotypes[1], 1);
        }
      } else {
        fail("Can only test on 1 or 2 haplotypes");
      }

      GeneReport geneReport = getContext().getGeneReport(gene);
      assertTrue(geneReport.isReportable(), "Not reportable: " + geneReport.getReporterDiplotypes());
      assertTrue(geneReport.getReporterDiplotypes().stream()
          .anyMatch(d -> d.makeLookupMap().equals(lookup)), "Lookup key " + lookup + " not found in lookup " + geneReport.getReporterDiplotypes().stream().map(Diplotype::makeLookupMap).toList());
    }

    private void testLookupByActivity(String gene, String activityScore) {
      GeneReport geneReport = getContext().getGeneReport(gene);
      assertTrue(geneReport.isReportable());
      assertTrue(geneReport.getReporterDiplotypes().stream()
          .allMatch(d -> d.printLookupKeys().equals(activityScore)));
    }

    /**
     * Check to see if all the given genes have been called by the matcher
     */
    private void testCalledByMatcher(String... genes) {
      assertTrue(genes != null && genes.length > 0);

      Arrays.stream(genes)
          .forEach(g -> assertTrue(getContext().getGeneReport(g).isCalled(), g + " is not called"));
    }

    private void testReportable(String... genes) {
      assertTrue(genes != null && genes.length > 0);
      Arrays.stream(genes)
          .forEach(g -> assertTrue(getContext().getGeneReport(g).isReportable(), g + " is not reportable"));
    }

    /**
     * Check to see if none of the given genes have been called by the matcher
     */
    private void testNotCalledByMatcher(String... genes) {
      assertTrue(genes != null && genes.length > 0);

      Arrays.stream(genes)
          .forEach(g -> assertFalse(getContext().getGeneReport(g).isCalled(), g + " is called"));
    }

    /**
     * Check to see if there is a matching recommendation for the given drug name
     * @param drugName a drug name that has recommendations
     * @param expectedCount the number of matching recommendations you expect
     */
    private void testMatchedGroups(String drugName, int expectedCount) {
      DrugReport drugReport = getContext().getDrugReport(drugName);
      assertEquals(expectedCount, drugReport.getMatchedGroupCount(),
          drugName + " does not have matching recommendation count of " + expectedCount + " (found " +
              drugReport.getMatchedGroupCount() + ")");
    }

    private void testMatchedGroups(String drugName, DataSource source, int expectedCount) {
      DrugReport drugReport = getContext().getDrugReport(drugName);
      assertEquals(
          expectedCount,
          drugReport.getGuidelines().stream().filter(g -> g.getSource() == source).mapToLong(g -> g.getAnnotationGroups().size()).sum(),
          drugName + " does not have matching recommendation count of " + expectedCount + " (found " +
              drugReport.getMatchedGroupCount() + ")");
    }

    private void testAnyMatchFromSource(String drugName, DataSource source) {
      DrugReport drugReport = getContext().getDrugReport(drugName);
      assertTrue(drugReport.getGuidelines().stream().anyMatch((g) -> g.getSource() == source && g.isMatched()),
          drugName + " does not have matching recommendation from " + source);
    }

    private void testNoMatchFromSource(String drugName, DataSource source) {
      DrugReport drugReport = getContext().getDrugReport(drugName);
      assertTrue(drugReport.getGuidelines().stream().noneMatch(r -> r.getSource() == source && r.isMatched()),
          drugName + " has a matching recommendation from " + source + " and expected none");
    }

    private void testMessageCountForDrug(String drugName, int messageCount) {
      DrugReport guideline = getContext().getDrugReport(drugName);
      assertEquals(messageCount, guideline.getMessages().size(),
          drugName + " expected " + messageCount + " messages and got " + guideline.getMessages());
    }
  }
}

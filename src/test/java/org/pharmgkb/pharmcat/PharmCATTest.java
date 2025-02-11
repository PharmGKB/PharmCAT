package org.pharmgkb.pharmcat;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Collection;
import java.util.Optional;
import java.util.function.Consumer;
import java.util.stream.Stream;
import org.checkerframework.checker.nullness.qual.Nullable;
import org.jsoup.Jsoup;
import org.jsoup.nodes.Document;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.TestInfo;
import org.pharmgkb.common.util.PathUtils;
import org.pharmgkb.pharmcat.haplotype.ResultSerializer;
import org.pharmgkb.pharmcat.haplotype.model.GeneCall;
import org.pharmgkb.pharmcat.haplotype.model.Result;
import org.pharmgkb.pharmcat.phenotype.Phenotyper;
import org.pharmgkb.pharmcat.reporter.MessageHelper;
import org.pharmgkb.pharmcat.reporter.model.DataSource;
import org.pharmgkb.pharmcat.reporter.model.result.CallSource;
import org.pharmgkb.pharmcat.reporter.model.result.GeneReport;
import org.pharmgkb.pharmcat.reporter.model.result.Haplotype;

import static com.github.stefanbirkner.systemlambda.SystemLambda.tapSystemOut;
import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.containsString;
import static org.junit.jupiter.api.Assertions.*;


/**
 * This is a JUnit test for {@link PharmCAT}.
 * This tests the CLI.
 *
 * @author Ryan Whaley
 */
class PharmCATTest {

  @AfterEach
  void deleteDirectory(TestInfo testInfo) {
    TestUtils.deleteTestOutputDirectory(testInfo);
  }


  @Test
  void basic() throws Exception {
    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/reference.vcf");
    Path refMatcherOutput = vcfFile.getParent().resolve("reference.match.json");
    Path refPhenotyperOutput = vcfFile.getParent().resolve("reference.phenotype.json");
    Path refReporterOutput = vcfFile.getParent().resolve("reference.report.html");

    // require VCF
    String systemOut = tapSystemOut(() -> PharmCAT.main(null));
    //System.out.println(systemOut);
    assertTrue(systemOut.contains("No input"));
    assertTrue(systemOut.contains("-vcf"));


    // standard full run - should output to the same directory as the VCF file
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
      assertNotNull(document.getElementById("ABCG2"));
      assertNull(document.getElementById("CYP2D6"));
      assertNotNull(document.getElementById("desflurane"));
      assertNull(document.getElementById("aripiprazole"));

    } finally {
      // always delete
      boolean save = TestUtils.isSaveTestOutput();
      TestUtils.setSaveTestOutput(false);
      TestUtils.deleteTestFiles(refMatcherOutput, refPhenotyperOutput, refReporterOutput);
      TestUtils.setSaveTestOutput(save);
    }
  }


  @Test
  void reference(TestInfo testInfo) throws Exception {
    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/reference.vcf");
    Path outputDir = TestUtils.getTestOutputDir(testInfo, true);
    doReference(testInfo, new String[]{
        "-vcf", vcfFile.toString(),
        "-o", outputDir.toString(),
    }, document -> {
      assertNotNull(document.getElementById("desflurane"));
      assertEquals(1, document.getElementsByClass("cpic-guideline-desflurane").size());
      assertEquals(1, document.getElementsByClass("fda-label-desflurane").size());
      assertNotNull(document.getElementById("atorvastatin"));
      assertEquals(1, document.getElementsByClass("cpic-guideline-atorvastatin").size());
      assertEquals(1, document.getElementsByClass("dpwg-guideline-atorvastatin").size());
      assertEquals(1, document.getElementsByClass("fda-assoc-atorvastatin").size());
    });
  }

  private void doReference(TestInfo testInfo, String[] args, Consumer<Document> docChecker) throws Exception{
    Path outputDir = TestUtils.getTestOutputDir(testInfo, true);
    Path refMatcherOutput = outputDir.resolve("reference.match.json");
    Path refPhenotyperOutput = outputDir.resolve("reference.phenotype.json");
    Path refReporterOutput = outputDir.resolve("reference.report.html");

    try {
      String systemOut = tapSystemOut(() -> PharmCAT.main(args));
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
      assertNotNull(document.getElementById("ABCG2"));
      assertNull(document.getElementById("CYP2D6"));
      docChecker.accept(document);

    } finally {
      TestUtils.deleteTestFiles(refMatcherOutput, refPhenotyperOutput, refReporterOutput);
    }
  }


  @Test
  void reference_cpicAndDpwgOnly(TestInfo testInfo) throws Exception {
    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/reference.vcf");
    Path outputDir = TestUtils.getTestOutputDir(testInfo, true);
    doReference(testInfo, new String[]{
        "-vcf", vcfFile.toString(),
        "-o", outputDir.toString(),
        "-rs", "CPIC,DPWG"
    }, document -> {
      assertNotNull(document.getElementById("desflurane"));
      assertEquals(1, document.getElementsByClass("cpic-guideline-desflurane").size());
      assertEquals(0, document.getElementsByClass("fda-label-desflurane").size());
      assertNotNull(document.getElementById("atorvastatin"));
      assertEquals(1, document.getElementsByClass("cpic-guideline-atorvastatin").size());
      assertEquals(1, document.getElementsByClass("dpwg-guideline-atorvastatin").size());
      assertEquals(0, document.getElementsByClass("fda-assoc-atorvastatin").size());
    });
  }


  @Test
  void reference_fdaOnly(TestInfo testInfo) throws Exception {
    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/reference.vcf");
    Path outputDir = TestUtils.getTestOutputDir(testInfo, true);
    doReference(testInfo, new String[]{
        "-vcf", vcfFile.toString(),
        "-o", outputDir.toString(),
        "-rs", "FDA"
    }, document -> {
      assertNotNull(document.getElementById("desflurane"));
      assertEquals(0, document.getElementsByClass("cpic-guideline-desflurane").size());
      assertEquals(1, document.getElementsByClass("fda-label-desflurane").size());
      assertNull(document.getElementById("atorvastatin"));
    });
  }

  @Test
  void referenceTsvOnly(TestInfo testInfo) throws Exception {
    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/reference.vcf");
    Path outputDir = TestUtils.getTestOutputDir(testInfo, true);

    Path refMatcherOutput = outputDir.resolve("reference.match.json");
    Path refPhenotyperOutput = outputDir.resolve("reference.phenotype.json");
    Path refReporterOutput = outputDir.resolve("reference.report.html");
    Path refReporterTsvOutput = outputDir.resolve("reference.report.tsv");

    try {
      String systemOut = tapSystemOut(() -> PharmCAT.main(new String[] {
          "-vcf", vcfFile.toString(),
          "-o", outputDir.toString(),
          "-del", "-reporterCallsOnlyTsv"
      }));
      System.out.println(systemOut);
      assertTrue(systemOut.contains("Done."));
      assertFalse(Files.exists(refMatcherOutput));
      assertFalse(Files.exists(refPhenotyperOutput));
      assertFalse(Files.exists(refReporterOutput));
      assertTrue(Files.exists(refReporterTsvOutput));

    } finally {
      TestUtils.deleteTestFiles(refMatcherOutput, refPhenotyperOutput, refReporterOutput);
    }
  }


  /**
   * An example run with CYP2D6 research mode enabled.
   * <p>NOTE: since research mode is enabled you will not get output from the reporter.</p>
   */
  @Test
  void callCyp2d6(TestInfo testInfo) throws Exception {
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
      assertTrue(Files.notExists(refReporterOutput));

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

    } finally {
      TestUtils.deleteTestFiles(outputDir);
    }
  }

  @Test
  void callCyp2d6WithOverlappingOutsideCall(TestInfo testInfo) throws Exception {
    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/reference.vcf");
    Path outsideCallFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/PharmCATTest-cyp2d6.tsv");

    Path outputDir = TestUtils.getTestOutputDir(testInfo, true);
    Path matcherOutput = outputDir.resolve("reference.match.json");
    Path phenotyperOutput = outputDir.resolve("reference.phenotype.json");
    Path reporterOutput = outputDir.resolve("reference.report.html");

    try {
      String systemOut = tapSystemOut(() -> PharmCAT.main(new String[]{
          "-vcf", vcfFile.toString(),
          "-po", outsideCallFile.toString(),
          "-o", outputDir.toString(),
          "-matcher", "-research", "cyp2d6"
      }));
      System.out.println(systemOut);
      assertTrue(Files.exists(matcherOutput));
      assertFalse(Files.exists(phenotyperOutput));
      assertFalse(Files.exists(reporterOutput));

    } finally {
      TestUtils.deleteTestFiles(outputDir);
    }
  }

  @Test
  void outsideCallsOnly(TestInfo testInfo) throws Exception {
    Path outsideCallFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/PharmCATTest-cyp2d6.tsv");

    Path outputDir = TestUtils.getTestOutputDir(testInfo, true);
    String baseFilename = TestUtils.getFullTestName(testInfo);
    Path matcherOutput = outputDir.resolve(baseFilename + ".match.json");
    Path phenotyperOutput = outputDir.resolve(baseFilename + ".phenotype.json");
    Path reporterOutput = outputDir.resolve(baseFilename + ".report.html");

    String systemOut = tapSystemOut(() -> PharmCAT.main(new String[] {
        "-o", outputDir.toString(),
        "-phenotyper",
        "-po", outsideCallFile.toString(),
        "-bf", baseFilename,
    }));
    System.out.println(systemOut);
    assertTrue(systemOut.contains("Done."));
    assertFalse(Files.exists(matcherOutput));
    assertTrue(Files.exists(phenotyperOutput));
    assertFalse(Files.exists(reporterOutput));

    validateCyp2d6OutsideCallOutput(phenotyperOutput);
  }

  @Test
  void outsideCallsNoRecs(TestInfo testInfo) throws Exception {
    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/reference.vcf");
    Path outsideCallFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/PharmCATTest-outsideCallsNoRecs.tsv");

    Path outputDir = TestUtils.getTestOutputDir(testInfo, true);
    Path matcherOutput = outputDir.resolve("reference.match.json");
    Path phenotyperOutput = outputDir.resolve("reference.phenotype.json");
    Path reporterOutput = outputDir.resolve("reference.report.html");

    try {
      String systemOut = tapSystemOut(() -> PharmCAT.main(new String[] {
          "-vcf", vcfFile.toString(),
          "-po", outsideCallFile.toString(),
          "-o", outputDir.toString(),
      }));
      System.out.println(systemOut);
      assertTrue(systemOut.contains("Done."));
      assertTrue(Files.exists(matcherOutput));
      assertTrue(Files.exists(phenotyperOutput));
      assertTrue(Files.exists(reporterOutput));

      Collection<GeneReport> reports = Phenotyper.read(phenotyperOutput).getGeneReports().get(DataSource.CPIC)
          .values();
      Optional<GeneReport> grOpt = reports.stream()
          .filter(gr -> gr.getGene().equals("IFNL3"))
          .findFirst();
      assertTrue(grOpt.isPresent());
      assertTrue(grOpt.get().isCalled());
      assertTrue(grOpt.get().isOutsideCall());

      Document document = Jsoup.parse(reporterOutput.toFile());
      assertEquals(1,
          document.select(".gene.IFNL3_4 .alert-warning." + MessageHelper.MSG_OUTSIDE_CALL).size());

    } finally {
      TestUtils.deleteTestFiles(outputDir);
    }
  }


  @Test
  void outsideCallsSuballeles(TestInfo testInfo) throws Exception {
    Path outsideCallFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/PharmCATTest-outsideCallsSuballeles.tsv");

    Path outputDir = TestUtils.getTestOutputDir(testInfo, true);
    Path matcherOutput = outputDir.resolve("PharmCATTest-outsideCallsSuballeles.match.json");
    Path phenotyperOutput = outputDir.resolve("PharmCATTest-outsideCallsSuballeles.phenotype.json");
    Path reporterOutput = outputDir.resolve("PharmCATTest-outsideCallsSuballeles.report.html");

    try {
      String systemOut = tapSystemOut(() -> PharmCAT.main(new String[] {
          "-o", outputDir.toString(),
          "-phenotyper",
          "-po", outsideCallFile.toString()
      }));
      System.out.println(systemOut);
      assertTrue(systemOut.contains("Done."));
      assertTrue(systemOut.contains("WARNING: PharmCAT does not support sub-alleles for CYP2D6."));
      assertTrue(systemOut.contains("WARNING: PharmCAT does not support sub-alleles for HLA-A."));
      assertTrue(systemOut.contains("WARNING: Converting outside call for HLA-B from 'B*04:05', to '*04:05'."));

      assertFalse(Files.exists(matcherOutput));
      assertTrue(Files.exists(phenotyperOutput));
      assertFalse(Files.exists(reporterOutput));

      validateOutsideCallOutput(phenotyperOutput, "CYP2D6", "*1/*2");
      validateOutsideCallOutput(phenotyperOutput, "HLA-A", "*02:02/*03:02");
      validateOutsideCallOutput(phenotyperOutput, "HLA-B", "*02:01/*04:05");
    } finally {
      TestUtils.deleteTestFiles(outputDir);
    }
  }



  /**
   * This test is dependent on CYP2C19 definition.
   * It may fail if CYP2C19 definition changes too much.
   */
  @Test
  void matchAllFlag(TestInfo testInfo) throws Exception {
    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/PharmCATTest-cyp2c19MissingPositions.vcf");

    Path outputDir = TestUtils.getTestOutputDir(testInfo, true);
    String baseFilename = TestUtils.getTestName(testInfo);
    Path matcherOutput = outputDir.resolve(baseFilename + BaseConfig.MATCHER_SUFFIX + ".json");
    Path phenotyperOutput = outputDir.resolve(baseFilename + BaseConfig.PHENOTYPER_SUFFIX + ".json");
    Path reporterOutput = outputDir.resolve(baseFilename + BaseConfig.REPORTER_SUFFIX + ".html");

    // matcher only, expecting 1 CYP2C19 matches
    try {
      String systemOut = tapSystemOut(() -> PharmCAT.main(new String[] {
          "-vcf", vcfFile.toString(),
          "-matcher",
          "-o", outputDir.toString(),
          "-bf", baseFilename,
      }));
      //System.out.println(systemOut);
      assertThat(systemOut, containsString("Done"));
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
    assertTrue(gc.getDiplotypes().size() > 50, "Expecting more than 50, found " + gc.getDiplotypes().size());
  }


  /**
   * This test ensures that the output of the phenotyper is consistent regardless of whether you run the matcher first
   * in one process then the phenotyper in a second process versus if you run the matcher and phenotyper together in
   * one process.
   */
  @Test
  void consistentOutput(TestInfo testInfo) throws Exception {
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

      StringBuilder singlePhenoJson = new StringBuilder();
      try (Stream<String> lines = Files.lines(singlesPhenotyperOutput)) {
        // Process each line
        lines.filter(l -> !l.contains("\"timestamp\":"))
            .forEach(l -> {
              singlePhenoJson.append(l)
                  .append("\n");
            });
      }
      StringBuilder doublePhenoJson = new StringBuilder();
      try (Stream<String> lines = Files.lines(doublePhenotyperOutput)) {
        lines.filter(l -> !l.contains("\"timestamp\":"))
            .forEach(l -> {
              doublePhenoJson.append(l)
                  .append("\n");
            });
      }

      assertEquals(singlePhenoJson.toString(), doublePhenoJson.toString());
    } finally {
      TestUtils.deleteTestFiles(outputDir);
    }
  }


  @Test
  void multisample(TestInfo testInfo) throws Exception {
    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/VcfSampleReaderTest.vcf");

    Path outputDir = TestUtils.getTestOutputDir(testInfo, true);
    Path matcherOutput1 = outputDir.resolve("VcfSampleReaderTest.Sample_1.match.json");
    Path phenotyperOutput1 = outputDir.resolve("VcfSampleReaderTest.Sample_1.phenotype.json");
    Path reporterOutput1 = outputDir.resolve("VcfSampleReaderTest.Sample_1.report.html");
    Path matcherOutput2 = outputDir.resolve("VcfSampleReaderTest.Sample_2.match.json");
    Path phenotyperOutput2 = outputDir.resolve("VcfSampleReaderTest.Sample_2.phenotype.json");
    Path reporterOutput2 = outputDir.resolve("VcfSampleReaderTest.Sample_2.report.html");

    try {
      String systemOut = tapSystemOut(() -> PharmCAT.main(new String[] {
          "-vcf", vcfFile.toString(),
          "-o", outputDir.toString(),
      }));
      System.out.println(systemOut);
      assertTrue(systemOut.contains("Done."));

      assertTrue(Files.exists(matcherOutput1));
      assertTrue(Files.exists(phenotyperOutput1));
      assertTrue(Files.exists(reporterOutput1));

      assertTrue(Files.exists(matcherOutput2));
      assertTrue(Files.exists(phenotyperOutput2));
      assertTrue(Files.exists(reporterOutput2));

    } finally {
      TestUtils.deleteTestFiles(outputDir);
    }
  }

  @Test
  void multipleOutsideCallFiles(TestInfo testInfo) throws Exception {
    Path outsideCallFile1 = PathUtils.getPathToResource("org/pharmgkb/pharmcat/PharmCATTest-cyp2d6.tsv");
    Path outsideCallFile2 = PathUtils.getPathToResource("org/pharmgkb/pharmcat/PharmCATTest-outsideCallsNoRecs.tsv");
    Path outputDir = TestUtils.getTestOutputDir(testInfo, true);

    String systemOut = tapSystemOut(() -> PharmCAT.main(new String[] {
        "-phenotyper",
        "-reporter",
        "-reporterHtml",
        "-reporterJson",
        "-po", outsideCallFile1.toString(),
        "-po", outsideCallFile2.toString(),
        "-o", outputDir.toString(),
    }));
    assertTrue(systemOut.contains("Done."));

    // file names should be based on the first outside call file
    assertTrue(Files.exists(outputDir.resolve("PharmCATTest-cyp2d6.phenotype.json")));
    assertTrue(Files.exists(outputDir.resolve("PharmCATTest-cyp2d6.report.json")));
    assertTrue(Files.exists(outputDir.resolve("PharmCATTest-cyp2d6.report.html")));

    // file names should NOT be based on the second outside call file
    assertFalse(Files.exists(outputDir.resolve("PharmCATTest-outsideCallsNoRecs.phenotype.json")));
    assertFalse(Files.exists(outputDir.resolve("PharmCATTest-outsideCallsNoRecs.report.json")));
    assertFalse(Files.exists(outputDir.resolve("PharmCATTest-outsideCallsNoRecs.report.html")));

    Phenotyper phenotyper = Phenotyper.read(outputDir.resolve("PharmCATTest-cyp2d6.phenotype.json"));
    checkOutsideDiplotype(phenotyper.findGeneReport(DataSource.CPIC, "CYP2D6").orElse(null),
        "*3", "*4");
    checkOutsideDiplotype(phenotyper.findGeneReport(DataSource.CPIC, "CYP4F2").orElse(null),
        "*1", "*3");
    checkOutsideDiplotype(phenotyper.findGeneReport(DataSource.CPIC, "IFNL3").orElse(null),
        "rs12979860 variant (T)", "rs12979860 variant (T)");
  }

  public static void checkOutsideDiplotype(@Nullable GeneReport report, String allele1, String allele2) {
    assertNotNull(report);
    assertEquals(CallSource.OUTSIDE, report.getCallSource());
    Haplotype haplotype = report.getSourceDiplotypes().first().getAllele1();
    assertNotNull(haplotype);
    assertEquals(allele1, haplotype.getName());
    haplotype = report.getSourceDiplotypes().first().getAllele2();
    assertNotNull(haplotype);
    assertEquals(allele2, haplotype.getName());
  }


  public static void validateCyp2d6OutsideCallOutput(Path phenotyperOutput) throws IOException {
    validateOutsideCallOutput(phenotyperOutput, "CYP2D6", null);;
  }


  public static void validateOutsideCallOutput(Path phenotyperOutput, String gene, @Nullable String diplotype)
      throws IOException {
    Collection<GeneReport> reports = Phenotyper.read(phenotyperOutput).getGeneReports().get(DataSource.CPIC)
        .values();
    Optional<GeneReport> grOpt = reports.stream()
        .filter(gr -> gr.getGene().equals(gene))
        .findFirst();
    assertTrue(grOpt.isPresent());
    assertTrue(grOpt.get().isCalled());
    assertTrue(grOpt.get().isOutsideCall());
    if (diplotype != null) {
      assertEquals(diplotype, grOpt.get().getSourceDiplotypes().first().getLabel());
    }
  }
}

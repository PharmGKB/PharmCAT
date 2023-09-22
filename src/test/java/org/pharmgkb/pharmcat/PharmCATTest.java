package org.pharmgkb.pharmcat;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Collection;
import java.util.Optional;
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
import org.pharmgkb.pharmcat.reporter.model.DataSource;
import org.pharmgkb.pharmcat.reporter.model.result.GeneReport;

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
      assertNotNull(document.getElementById("ABCG2"));
      assertNull(document.getElementById("CYP2D6"));
      assertNotNull(document.getElementById("desflurane"));

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

      validateCyp2d6OutsideCallOutput(phenotyperOutput);

    } finally {
      TestUtils.deleteTestFiles(outputDir);
    }
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
      assertEquals(1, document.select(".gene.IFNL3_4 .alert-warning.pcat-outside-call").size());

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

      assertEquals(Files.readString(singlesPhenotyperOutput), Files.readString(doublePhenotyperOutput));
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


  public static void validateCyp2d6OutsideCallOutput(Path phenotyperOutput) throws IOException {
    Collection<GeneReport> reports = Phenotyper.read(phenotyperOutput).getGeneReports().get(DataSource.CPIC)
        .values();
    Optional<GeneReport> grOpt = reports.stream()
        .filter(gr -> gr.getGene().equals("CYP2D6"))
        .findFirst();
    assertTrue(grOpt.isPresent());
    assertTrue(grOpt.get().isCalled());
    assertTrue(grOpt.get().isOutsideCall());
  }
}

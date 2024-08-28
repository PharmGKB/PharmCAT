package org.pharmgkb.pharmcat;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Collection;
import java.util.Optional;
import java.util.function.Consumer;
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
import org.pharmgkb.pharmcat.reporter.model.DataSource;
import org.pharmgkb.pharmcat.reporter.model.result.GeneReport;

import static com.github.stefanbirkner.systemlambda.SystemLambda.tapSystemOut;
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

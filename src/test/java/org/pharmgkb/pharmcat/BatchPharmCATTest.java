package org.pharmgkb.pharmcat;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import com.google.common.base.Stopwatch;
import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;
import org.apache.commons.io.FilenameUtils;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.TestInfo;
import org.pharmgkb.common.util.PathUtils;
import org.pharmgkb.pharmcat.haplotype.VcfSampleReader;
import org.pharmgkb.pharmcat.phenotype.Phenotyper;
import org.pharmgkb.pharmcat.reporter.model.DataSource;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.*;
import static org.junit.jupiter.api.Assertions.*;
import static uk.org.webcompere.systemstubs.SystemStubs.tapSystemErr;
import static uk.org.webcompere.systemstubs.SystemStubs.tapSystemOut;


/**
 * This is a JUnit test for {@link BatchPharmCAT}.
 * This tests the CLI.
 *
 * @author Mark Woon
 */
class BatchPharmCATTest {

  @AfterEach
  void deleteDirectory(TestInfo testInfo) {
    TestUtils.deleteTestOutputDirectory(testInfo);
  }


  @Test
  void noArgs() throws Exception {
    // require input directory
    String systemErr = tapSystemErr(() -> BatchPharmCAT.main(null));
    //System.out.println(systemErr);
    assertThat(systemErr, containsString("Missing input"));
  }

  @Test
  void noFiles(TestInfo testInfo) throws Exception {
    Path tmpDir = TestUtils.getTestOutputDir(testInfo, true);
    String systemOut = tapSystemOut(() -> BatchPharmCAT.main(new String[] {
        "-i", tmpDir.toString(),
    }));
    System.out.println(systemOut);
    assertThat(systemOut, containsString("No input files"));
  }

  @Test
  void simpleRun(TestInfo testInfo) throws Exception {
    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/reference.vcf");

    Path tmpDir = TestUtils.getTestOutputDir(testInfo, true);
    copyFiles(tmpDir, vcfFile);

    String systemOut = tapSystemOut(() -> BatchPharmCAT.main(new String[] {
        "-i", tmpDir.toString(),
        "-cp", "3"
    }));
    System.out.println(systemOut);
    assertThat(systemOut, containsString("Done."));
    assertThat(systemOut, not(containsString("FAIL")));
    // max processes value is capped to the number of samples, so capped to 1, which is not shown
    assertThat(systemOut, not(matchesPattern("maximum of \\d+ processes")));
    checkForOutputFiles(tmpDir, vcfFile);
  }


  @Test
  void compressed(TestInfo testInfo) throws Exception {
    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/multisample.vcf.bgz");

    Path tmpDir = TestUtils.getTestOutputDir(testInfo, true);
    copyFiles(tmpDir, vcfFile);

    String systemOut = tapSystemOut(() -> BatchPharmCAT.main(new String[] {
        "-vcf", tmpDir.resolve(vcfFile.getFileName()).toString(),
    }));
    System.out.println(systemOut);
    assertThat(systemOut, containsString("Done."));
    assertThat(systemOut, not(containsString("FAIL")));
    if (!TestUtils.isContinuousIntegration()) {
      // max processes is capped to number of samples
      // don't test this in CI because no way of guaranteeing # of processors
      assertThat(systemOut, containsString("maximum of 2 processes"));
    }

    checkForOutputFiles(tmpDir, vcfFile);

    Path w1 = tmpDir.resolve("multisample.Sample_1.match_warnings.txt");
    assertTrue(Files.exists(w1), "Missing " + w1);
    Path w2 = tmpDir.resolve("multisample.Sample_2.match_warnings.txt");
    assertTrue(Files.exists(w1), "Missing " + w2);
  }


  @Test
  void twoSamples(TestInfo testInfo) throws Exception {
    Path[] vcfFiles = new Path[] {
        PathUtils.getPathToResource("org/pharmgkb/pharmcat/Sample_1.preprocessed.vcf"),
        PathUtils.getPathToResource("org/pharmgkb/pharmcat/Sample_2.preprocessed.vcf"),
    };

    Path tmpDir = TestUtils.getTestOutputDir(testInfo, true);
    copyFiles(tmpDir, vcfFiles);

    String systemOut = tapSystemOut(() -> BatchPharmCAT.main(new String[] {
        "-i", tmpDir.toString(),
    }));
    System.out.println(systemOut);
    assertThat(systemOut, containsString("Queueing up 2 samples"));
    assertThat(systemOut, containsString("Done."));
    assertThat(systemOut, not(containsString("FAIL")));
    checkForOutputFiles(tmpDir, vcfFiles);
  }

  @Test
  void twoOutsideCallsForOneSample(TestInfo testInfo) throws Exception {
    Path[] vcfFiles = new Path[] {
        PathUtils.getPathToResource("org/pharmgkb/pharmcat/Sample_1.preprocessed.vcf"),
        PathUtils.getPathToResource("org/pharmgkb/pharmcat/Sample_2.preprocessed.vcf"),
    };

    Path tmpDir = TestUtils.getTestOutputDir(testInfo, true);
    copyFiles(tmpDir, vcfFiles);

    Path outsideFileS2O1 = tmpDir.resolve("Sample_2.outside1.tsv");
    Files.copy(PathUtils.getPathToResource("org/pharmgkb/pharmcat/PharmCATTest-cyp2d6.tsv"), outsideFileS2O1);

    Path outsideFileS2O2 = tmpDir.resolve("Sample_2.outside2.tsv");
    Files.copy(PathUtils.getPathToResource("org/pharmgkb/pharmcat/PharmCATTest-outsideCallsNoRecs.tsv"), outsideFileS2O2);

    String systemOut = tapSystemOut(() -> BatchPharmCAT.main(new String[] {
        "-i", tmpDir.toString(),
    }));
    System.out.println(systemOut);
    assertThat(systemOut, containsString("Found 2 VCF files"));
    // this should always come out as 2 samples queued, the two outside call files are attached to one sample
    assertThat(systemOut, containsString("Queueing up 2 samples"));
    assertThat(systemOut, containsString("Done."));
    assertThat(systemOut, not(containsString("FAIL")));
    checkForOutputFiles(tmpDir, vcfFiles);

    Phenotyper phenotyper = Phenotyper.read(tmpDir.resolve("Sample_2.phenotype.json"));
    PharmCATTest.checkOutsideDiplotype(phenotyper.findGeneReport(DataSource.CPIC, "CYP2D6").orElse(null),
        "*3", "*4");
    PharmCATTest.checkOutsideDiplotype(phenotyper.findGeneReport(DataSource.CPIC, "CYP4F2").orElse(null),
        "*1", "*3");
    PharmCATTest.checkOutsideDiplotype(phenotyper.findGeneReport(DataSource.CPIC, "IFNL3").orElse(null),
        "rs12979860 variant (T)", "rs12979860 variant (T)");
  }


  @Test
  void sixSamples(TestInfo testInfo) throws Exception {
    Path[] vcfFiles = new Path[] {
        PathUtils.getPathToResource("org/pharmgkb/pharmcat/Sample_1.preprocessed.vcf"),
        PathUtils.getPathToResource("org/pharmgkb/pharmcat/Sample_2.preprocessed.vcf"),
        PathUtils.getPathToResource("org/pharmgkb/pharmcat/PharmCATTest-cyp2c19MissingPositions.vcf"),
    };

    Path tmpDir = TestUtils.getTestOutputDir(testInfo, true);
    copyFiles(tmpDir, vcfFiles);

    Path outsideFile1 = tmpDir.resolve("Sample_1.outside.tsv");
    Path matchFile3 = tmpDir.resolve("Sample_3.match.json");
    Path outsideFile4 = tmpDir.resolve("Sample_4.outside.tsv");
    Path phenotypeFile5 = tmpDir.resolve("Sample5.phenotype.json");
    Files.copy(PathUtils.getPathToResource("org/pharmgkb/pharmcat/PharmCATTest-cyp2d6.tsv"), outsideFile1);
    Files.copy(PathUtils.getPathToResource("org/pharmgkb/pharmcat/PharmCATTest-cyp2d6.tsv"), outsideFile4);
    Files.copy(PathUtils.getPathToResource("org/pharmgkb/pharmcat/Sample_1.match.json"), matchFile3);
    Files.copy(PathUtils.getPathToResource("org/pharmgkb/pharmcat/Sample_1.phenotype.json"), phenotypeFile5);

    String systemOut = tapSystemOut(() -> BatchPharmCAT.main(new String[] {
        "-i", tmpDir.toString(),
    }));
    System.out.println(systemOut);
    assertThat(systemOut, containsString("Done."));
    assertThat(systemOut, not(containsString("FAIL")));
    assertThat(systemOut, containsString("Found 3 VCF files"));
    assertThat(systemOut, containsString("Found 1 independent phenotyper input file"));
    assertThat(systemOut, containsString("Found 1 independent phenotyper outside call file"));
    assertThat(systemOut, containsString("Warning: lone outside call file"));
    assertThat(systemOut, containsString("Found 1 independent reporter input file"));
    assertThat(systemOut, containsString("Queueing up 6 samples"));

    List<Path> allInputs = new ArrayList<>(Arrays.asList(vcfFiles));
    allInputs.add(outsideFile1);
    allInputs.add(matchFile3);
    allInputs.add(outsideFile4);
    allInputs.add(phenotypeFile5);
    checkForOutputFiles(tmpDir, allInputs.toArray(new Path[0]));

    PharmCATTest.validateCyp2d6OutsideCallOutput(tmpDir.resolve("Sample_1.phenotype.json"));
    PharmCATTest.validateCyp2d6OutsideCallOutput(tmpDir.resolve("Sample_4.phenotype.json"));
  }


  @Test
  void multisample(TestInfo testInfo) throws Exception {
    Path na18526Vcf  = PathUtils.getPathToResource("org/pharmgkb/pharmcat/PharmCATTest-cyp2c19MissingPositions.vcf");
    Path multisampleVcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/VcfSampleReaderTest-multisample.vcf");
    Path[] vcfFiles = new Path[] {
        multisampleVcfFile,
        na18526Vcf
    };

    Path tmpDir = TestUtils.getTestOutputDir(testInfo, true);
    copyFiles(tmpDir, vcfFiles);

    Path outsideFile1 = tmpDir.resolve("Sample_1.outside.tsv");
    Path matchFile3 = tmpDir.resolve("Sample_3.match.json");
    Path outsideFile4 = tmpDir.resolve("Sample_4.outside.tsv");
    Path phenotypeFile5 = tmpDir.resolve("Sample5.phenotype.json");
    Files.copy(PathUtils.getPathToResource("org/pharmgkb/pharmcat/PharmCATTest-cyp2d6.tsv"), outsideFile1);
    Files.copy(PathUtils.getPathToResource("org/pharmgkb/pharmcat/PharmCATTest-cyp2d6.tsv"), outsideFile4);
    Files.copy(PathUtils.getPathToResource("org/pharmgkb/pharmcat/Sample_1.match.json"), matchFile3);
    Files.copy(PathUtils.getPathToResource("org/pharmgkb/pharmcat/Sample_1.phenotype.json"), phenotypeFile5);

    String systemOut = tapSystemOut(() -> BatchPharmCAT.main(new String[] {
        "-i", tmpDir.toString(),
        "-v"
    }));
    System.out.println(systemOut);
    assertThat(systemOut, containsString("Done."));
    assertThat(systemOut, not(containsString("FAIL")));
    assertThat(systemOut, containsString("Found 2 VCF files"));
    assertThat(systemOut, containsString("Found 1 independent phenotyper input file"));
    assertThat(systemOut, containsString("Found 1 independent phenotyper outside call file"));
    assertThat(systemOut, containsString("Warning: lone outside call file"));
    assertThat(systemOut, containsString("Found 1 independent reporter input file"));
    assertThat(systemOut, containsString("Queueing up 6 samples"));
    if (!TestUtils.isContinuousIntegration()) {
      // max processes is capped to number of samples
      // don't test this in CI because no way of guaranteeing # of processors
      assertThat(systemOut, containsString("maximum of 6 processes"));
    }

    List<Path> allInputs = new ArrayList<>();
    allInputs.add(multisampleVcfFile);
    allInputs.add(na18526Vcf);
    // don't add outside call because multisample VCFs handle base names differently
    //allInputs.add(outsideFile1);
    allInputs.add(matchFile3);
    allInputs.add(outsideFile4);
    allInputs.add(phenotypeFile5);
    checkForOutputFiles(tmpDir, allInputs.toArray(new Path[0]));

    PharmCATTest.validateCyp2d6OutsideCallOutput(tmpDir.resolve("VcfSampleReaderTest-multisample.Sample_1.phenotype.json"));
    PharmCATTest.validateCyp2d6OutsideCallOutput(tmpDir.resolve("Sample_4.phenotype.json"));
  }


  @Test
  void multisampleRestricted(TestInfo testInfo) throws Exception {
    Path na18526Vcf  = PathUtils.getPathToResource("org/pharmgkb/pharmcat/PharmCATTest-cyp2c19MissingPositions.vcf");
    Path multisampleVcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/VcfSampleReaderTest-multisample.vcf");
    Path[] vcfFiles = new Path[] {
        multisampleVcfFile,
        na18526Vcf
    };

    Path tmpDir = TestUtils.getTestOutputDir(testInfo, true);
    copyFiles(tmpDir, vcfFiles);

    Path outsideFile1 = tmpDir.resolve("Sample_1.outside.tsv");
    Path matchFile3 = tmpDir.resolve("Sample_3.match.json");
    Path outsideFile4 = tmpDir.resolve("Sample_4.outside.tsv");
    Path phenotypeFile5 = tmpDir.resolve("Sample5.phenotype.json");
    Files.copy(PathUtils.getPathToResource("org/pharmgkb/pharmcat/PharmCATTest-cyp2d6.tsv"), outsideFile1);
    Files.copy(PathUtils.getPathToResource("org/pharmgkb/pharmcat/PharmCATTest-cyp2d6.tsv"), outsideFile4);
    Files.copy(PathUtils.getPathToResource("org/pharmgkb/pharmcat/Sample_1.match.json"), matchFile3);
    Files.copy(PathUtils.getPathToResource("org/pharmgkb/pharmcat/Sample_1.phenotype.json"), phenotypeFile5);

    String systemOut = tapSystemOut(() -> BatchPharmCAT.main(new String[] {
        "-i", tmpDir.toString(),
        "-cp", "3"
    }));
    System.out.println(systemOut);
    assertThat(systemOut, containsString("Done."));
    assertThat(systemOut, not(containsString("FAIL")));
    assertThat(systemOut, containsString("Found 2 VCF files"));
    assertThat(systemOut, containsString("Found 1 independent phenotyper input file"));
    assertThat(systemOut, containsString("Found 1 independent phenotyper outside call file"));
    assertThat(systemOut, containsString("Warning: lone outside call file"));
    assertThat(systemOut, containsString("Found 1 independent reporter input file"));
    assertThat(systemOut, containsString("Queueing up 6 samples"));
    if (!TestUtils.isContinuousIntegration()) {
      // max processes is lower than number of samples, so obey -cp
      // don't test this in CI because no way of guaranteeing # of processors
      assertThat(systemOut, containsString("maximum of 3 processes"));
    }

    List<Path> allInputs = new ArrayList<>();
    allInputs.add(multisampleVcfFile);
    allInputs.add(na18526Vcf);
    // don't add outside call because multisample VCFs handle base names differently
    //allInputs.add(outsideFile1);
    allInputs.add(matchFile3);
    allInputs.add(outsideFile4);
    allInputs.add(phenotypeFile5);
    checkForOutputFiles(tmpDir, allInputs.toArray(new Path[0]));

    PharmCATTest.validateCyp2d6OutsideCallOutput(tmpDir.resolve("VcfSampleReaderTest-multisample.Sample_1.phenotype.json"));
    PharmCATTest.validateCyp2d6OutsideCallOutput(tmpDir.resolve("Sample_4.phenotype.json"));
  }


  @Test
  void singleSampleIgnoreUnrelated(TestInfo testInfo) throws Exception {
    Path na18526Vcf  = PathUtils.getPathToResource("org/pharmgkb/pharmcat/PharmCATTest-cyp2c19MissingPositions.vcf");

    Path tmpDir = TestUtils.getTestOutputDir(testInfo, true);
    copyFiles(tmpDir, na18526Vcf);

    Path outsideFile1 = tmpDir.resolve("Sample_1.outside.tsv");
    Path matchFile3 = tmpDir.resolve("Sample_3.match.json");
    Path outsideFile4 = tmpDir.resolve("Sample_4.outside.tsv");
    Path phenotypeFile5 = tmpDir.resolve("Sample5.phenotype.json");
    Files.copy(PathUtils.getPathToResource("org/pharmgkb/pharmcat/PharmCATTest-cyp2d6.tsv"), outsideFile1);
    Files.copy(PathUtils.getPathToResource("org/pharmgkb/pharmcat/PharmCATTest-cyp2d6.tsv"), outsideFile4);
    Files.copy(PathUtils.getPathToResource("org/pharmgkb/pharmcat/Sample_1.match.json"), matchFile3);
    Files.copy(PathUtils.getPathToResource("org/pharmgkb/pharmcat/Sample_1.phenotype.json"), phenotypeFile5);

    String systemOut = tapSystemOut(() -> BatchPharmCAT.main(new String[] {
        "-vcf", na18526Vcf.toString(),
        "-o", tmpDir.toString(),
        "-cp", "3"
    }));
    System.out.println(systemOut);
    assertThat(systemOut, containsString("Done."));
    assertThat(systemOut, not(containsString("FAIL")));
    assertThat(systemOut, containsString("Found 1 VCF file"));
    assertThat(systemOut, not(containsString("independent phenotyper input file")));
    assertThat(systemOut, not(containsString("independent phenotyper outside call file")));
    assertThat(systemOut, not(containsString("lone outside call file")));
    assertThat(systemOut, not(containsString("independent reporter input file")));
    // max processes is higher than number of samples, should limit to 1, which is not shown
    assertThat(systemOut, not(matchesPattern("maximum of \\d+ processes")));

    checkForOutputFiles(tmpDir, na18526Vcf);
  }


  private void copyFiles(Path targetDir, Path... srcFiles) throws IOException {
    for (Path file : srcFiles) {
      Files.copy(file, targetDir.resolve(FilenameUtils.getName(file.toString())));
    }
  }


  private static final boolean sf_debugCheckOutput = false;

  private void checkForOutputFiles(Path dir, Path... inputFiles) throws IOException {
    Multimap<String, String> sampleMap = HashMultimap.create();
    for (Path file : inputFiles) {
      String baseName = BaseConfig.getBaseFilename(file);
      String extension = file.getFileName().toString().substring(baseName.length());
      if (sf_debugCheckOutput) {
        System.out.println("--> " + baseName + " -- " + extension);
      }
      boolean checked = false;
      if (extension.endsWith(".vcf") || extension.endsWith(".vcf.bgz") || extension.endsWith(".vcf.gz")) {
        if (sf_debugCheckOutput) {
          System.out.println("Checking .match.json");
        }
        VcfSampleReader reader = new VcfSampleReader(file);
        boolean singleSample = reader.getSamples().size() == 1;
        for (String sampleId : reader.getSamples()) {
          sampleMap.put(baseName, sampleId);
          Path f;
          if (singleSample) {
            f = dir.resolve(baseName + BaseConfig.MATCHER_SUFFIX + ".json");
          } else if (baseName.equals(sampleId)) {
            f = dir.resolve(baseName + BaseConfig.MATCHER_SUFFIX + ".json");
          } else {
            f = dir.resolve(baseName + "." + sampleId + BaseConfig.MATCHER_SUFFIX + ".json");
          }
          assertTrue(Files.exists(f), "Missing " + f);
        }
        checked = true;
      }
      if (extension.endsWith(".vcf") || extension.equals(BaseConfig.MATCHER_SUFFIX + ".json") ||
          BaseConfig.OUTSIDE_SUFFIX_PATTERN.matcher(extension).matches()) {
        if (sf_debugCheckOutput) {
          System.out.println("Checking .phenotype.json");
        }
        if (extension.endsWith(".vcf")) {
          boolean singleSample = sampleMap.get(baseName).size() == 1;
          for (String sampleId : sampleMap.get(baseName)) {
            Path f;
            if (singleSample) {
              f = dir.resolve(baseName + BaseConfig.MATCHER_SUFFIX + ".json");
            } else if (baseName.equals(sampleId)) {
              f = dir.resolve(sampleId + BaseConfig.PHENOTYPER_SUFFIX + ".json");
            } else {
              f = dir.resolve(baseName + "." + sampleId + BaseConfig.PHENOTYPER_SUFFIX + ".json");
            }
            assertTrue(Files.exists(f), "Missing " + f);
          }
        } else {
          Path f = dir.resolve(baseName + BaseConfig.PHENOTYPER_SUFFIX + ".json");
          assertTrue(Files.exists(f), "Missing " + f);
        }
        checked = true;
      }
      if (extension.endsWith(".vcf") || extension.equals(BaseConfig.MATCHER_SUFFIX + ".json") ||
          BaseConfig.OUTSIDE_EXTENSION_PATTERN.matcher(extension).matches() ||
          extension.equals(BaseConfig.PHENOTYPER_SUFFIX + ".json")) {
        if (sf_debugCheckOutput) {
          System.out.println("Checking .report.html");
        }
        if (extension.endsWith(".vcf")) {
          boolean singleSample = sampleMap.get(baseName).size() == 1;
          for (String sampleId : sampleMap.get(baseName)) {
            Path f;
            if (singleSample) {
              f = dir.resolve(baseName + BaseConfig.MATCHER_SUFFIX + ".json");
            } else if (baseName.equals(sampleId)) {
              f = dir.resolve(baseName + BaseConfig.REPORTER_SUFFIX + ".html");
            } else {
              f = dir.resolve(baseName + "." + sampleId + BaseConfig.REPORTER_SUFFIX + ".html");
            }
            assertTrue(Files.exists(f), "Missing " + f);
          }
        } else {
          Path f = dir.resolve(baseName + BaseConfig.REPORTER_SUFFIX + ".html");
          assertTrue(Files.exists(f), "Missing " + f);
        }
        checked = true;
      }
      if (!checked) {
        fail("Unrecognized extension: " + extension + " (" + file.getFileName() + ")");
      }
    }
  }

  @Test
  void compareToPharmcat100(TestInfo testInfo) throws Exception {
    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/100samples.vcf.bgz");

    Path batchTmpDir = TestUtils.getTestOutputDir(testInfo, "batch", true);
    copyFiles(batchTmpDir, vcfFile);

    String systemOut = tapSystemOut(() -> BatchPharmCAT.main(new String[] {
        "-vcf", batchTmpDir.resolve(vcfFile.getFileName().toString()).toString(),
        "-v"
    }));
    System.out.println(systemOut);
    assertThat(systemOut, containsString("Done."));
    assertThat(systemOut, not(containsString("FAIL")));


    System.out.println("-------------------------");

    Stopwatch stopwatch = Stopwatch.createStarted();
    Path serialTmpDir = TestUtils.getTestOutputDir(testInfo, "serial", true);
    copyFiles(serialTmpDir, vcfFile);
    systemOut = tapSystemOut(() -> PharmCAT.main(new String[] {
        "-vcf", serialTmpDir.resolve(vcfFile.getFileName().toString()).toString(),
        "-v"
    }));
    stopwatch.stop();
    System.out.println(systemOut);
    assertThat(systemOut, containsString("Done."));


    VcfSampleReader vsReader = new VcfSampleReader(vcfFile);
    List<String> samples = vsReader.getSamples();

    for (String sampleId : samples) {
      checkResult(batchTmpDir, serialTmpDir, "100samples", sampleId);
    }
  }

  private void checkResult(Path batchDir, Path serialDir, String basename, String sampleId) throws IOException {
    System.out.println("Checking " + sampleId);
    System.out.println("Checking matcher results");
    String batchMatchJson = readMatchJson(batchDir, basename, "match", sampleId);
    String serialMatchJson = readMatchJson(serialDir, basename, "match", sampleId);
    assertEquals(batchMatchJson, serialMatchJson, "Mismatch in " + sampleId);

    System.out.println("Checking phenotyper results");
    String batchPhenotypeJson = readMatchJson(batchDir, basename, "phenotype", sampleId);
    String serialPhenotypeJson = readMatchJson(serialDir, basename, "phenotype", sampleId);
    assertEquals(batchPhenotypeJson, serialPhenotypeJson, "Mismatch in " + sampleId);

    System.out.println("Checking reporter results");
    String batchReportHtml = Files.readString(batchDir.resolve("100samples." + sampleId + BaseConfig.REPORTER_SUFFIX + ".html"));
    String serialReportHtml = Files.readString(serialDir.resolve("100samples." + sampleId + BaseConfig.REPORTER_SUFFIX + ".html"));
    assertEquals(batchReportHtml, serialReportHtml, "Mismatch in " + sampleId);
  }


  private String readMatchJson(Path dir, String basename, String component, String sampleId) throws IOException {

    Path file = dir.resolve(basename + "." + sampleId + "." + component + ".json");
    assertTrue(Files.exists(file), "Cannot find " + file);
    StringBuilder builder = new StringBuilder();
    try (BufferedReader reader = Files.newBufferedReader(file)) {
      String line;
      while ((line = reader.readLine()) != null) {
        if (!line.contains("\"timestamp\":")) {
          builder.append(line)
              .append("\n");
        }
      }
    }
    return builder.toString();
  }
}

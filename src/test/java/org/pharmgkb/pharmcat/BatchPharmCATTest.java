package org.pharmgkb.pharmcat;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import org.apache.commons.io.FilenameUtils;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.TestInfo;
import org.pharmgkb.common.util.PathUtils;

import static com.github.stefanbirkner.systemlambda.SystemLambda.tapSystemErr;
import static com.github.stefanbirkner.systemlambda.SystemLambda.tapSystemOut;
import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.containsString;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.fail;


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
    // max processes is capped to number of samples
    assertThat(systemOut, containsString("maximum of 1 processes"));
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
    // max processes is capped to number of samples
    assertThat(systemOut, containsString("maximum of 2 processes"));

    checkForOutputFiles(tmpDir,
        tmpDir.resolve("multisample.Sample_1.vcf"),
        tmpDir.resolve("multisample.Sample_2.vcf"));

    Path w1 = tmpDir.resolve("multisample.Sample_1.matcher_warnings.txt");
    assertTrue(Files.exists(w1), "Missing " + w1);
    Path w2 = tmpDir.resolve("multisample.Sample_2.matcher_warnings.txt");
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
    checkForOutputFiles(tmpDir, vcfFiles);
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
    Path[] vcfFiles = new Path[] {
        PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/VcfSampleReaderTest.vcf"),
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
    allInputs.add(tmpDir.resolve("VcfSampleReaderTest.Sample_1.vcf"));
    allInputs.add(tmpDir.resolve("VcfSampleReaderTest.Sample_2.vcf"));
    allInputs.add(na18526Vcf);
    // don't add outside call because multisample VCFs handle base names differently
    //allInputs.add(outsideFile1);
    allInputs.add(matchFile3);
    allInputs.add(outsideFile4);
    allInputs.add(phenotypeFile5);
    checkForOutputFiles(tmpDir, allInputs.toArray(new Path[0]));

    PharmCATTest.validateCyp2d6OutsideCallOutput(tmpDir.resolve("VcfSampleReaderTest.Sample_1.phenotype.json"));
    PharmCATTest.validateCyp2d6OutsideCallOutput(tmpDir.resolve("Sample_4.phenotype.json"));
  }


  @Test
  void multisampleRestricted(TestInfo testInfo) throws Exception {
    Path na18526Vcf  = PathUtils.getPathToResource("org/pharmgkb/pharmcat/PharmCATTest-cyp2c19MissingPositions.vcf");
    Path[] vcfFiles = new Path[] {
        PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/VcfSampleReaderTest.vcf"),
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
    allInputs.add(tmpDir.resolve("VcfSampleReaderTest.Sample_1.vcf"));
    allInputs.add(tmpDir.resolve("VcfSampleReaderTest.Sample_2.vcf"));
    allInputs.add(na18526Vcf);
    // don't add outside call because multisample VCFs handle base names differently
    //allInputs.add(outsideFile1);
    allInputs.add(matchFile3);
    allInputs.add(outsideFile4);
    allInputs.add(phenotypeFile5);
    checkForOutputFiles(tmpDir, allInputs.toArray(new Path[0]));

    PharmCATTest.validateCyp2d6OutsideCallOutput(tmpDir.resolve("VcfSampleReaderTest.Sample_1.phenotype.json"));
    PharmCATTest.validateCyp2d6OutsideCallOutput(tmpDir.resolve("Sample_4.phenotype.json"));
  }


  private void copyFiles(Path targetDir, Path... srcFiles) throws IOException {
    for (Path file : srcFiles) {
      Files.copy(file, targetDir.resolve(FilenameUtils.getName(file.toString())));
    }
  }


  private static final boolean sf_debugCheckOutput = false;

  private void checkForOutputFiles(Path dir, Path... inputFiles) {
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
        Path f = dir.resolve(baseName + ".match.json");
        assertTrue(Files.exists(f), "Missing " + f);
        checked = true;
      }
      if (extension.endsWith(".vcf") || extension.equals(".match.json") || extension.equals(".outside.tsv")) {
        if (sf_debugCheckOutput) {
          System.out.println("Checking .phenotype.json");
        }
        Path f = dir.resolve(baseName + ".phenotype.json");
        assertTrue(Files.exists(f), "Missing " + f);
        checked = true;
      }
      if (extension.endsWith(".vcf") || extension.equals(".match.json") || extension.equals(".outside.tsv") ||
          extension.equals(".phenotype.json")) {
        if (sf_debugCheckOutput) {
          System.out.println("Checking .report.html");
        }
        Path f = dir.resolve(baseName + ".report.html");
        assertTrue(Files.exists(f), "Missing " + f);
        checked = true;
      }
      if (!checked) {
        fail("Unrecognized extension: " + extension + " (" + file.getFileName() + ")");
      }
    }
  }
}

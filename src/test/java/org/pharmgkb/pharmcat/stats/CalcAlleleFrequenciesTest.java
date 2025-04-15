package org.pharmgkb.pharmcat.stats;

import java.nio.file.Path;
import java.nio.file.Paths;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.TestInfo;
import org.pharmgkb.common.util.PathUtils;
import org.pharmgkb.pharmcat.PharmCAT;
import org.pharmgkb.pharmcat.ReportableException;
import org.pharmgkb.pharmcat.TestUtils;

import static org.junit.jupiter.api.Assertions.assertThrows;
import static uk.org.webcompere.systemstubs.SystemStubs.tapSystemOut;


/**
 * This is a JUnit test for {@link CalcAlleleFrequencies}.
 *
 * @author Mark Woon
 */
class CalcAlleleFrequenciesTest {

  @BeforeAll
  static void prepare() {
    //TestUtils.setSaveTestOutput(true);
  }

  @AfterEach
  void deleteDirectory(TestInfo testInfo) {
    TestUtils.deleteTestOutputDirectory(testInfo);
  }


  @Test
  void badInput() {

    Path outputDir = Paths.get("no/such/dir");
    assertThrows(ReportableException.class, () -> {
      String systemOut = tapSystemOut(() -> CalcAlleleFrequencies.main(new String[] {
          "-i", outputDir.toString()
      }));
      System.out.println(systemOut);
    });
  }


  @Test
  void test(TestInfo testInfo) throws Exception {
    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/VcfSampleReaderTest-multisample.vcf");
    Path saFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/VcfSampleReaderTest-multisample.sampleData.tsv");
    Path outputDir = TestUtils.getTestOutputDir(testInfo, true);

    PharmCAT.main(new String[] {
        "-vcf", vcfFile.toString(),
        "-o", outputDir.toString(),
        "-sm", saFile.toString(),
        "-reporterCallsOnlyTsv", "-del"
    });

    System.out.println();
    System.out.println("-----");
    System.out.println();

    String systemOut = tapSystemOut(() -> CalcAlleleFrequencies.main(new String[] {
        "-i", outputDir.toString()
    }));
    System.out.println(systemOut);

    System.out.println();
    System.out.println("-----");
    System.out.println();

    Path pivotOutputDir = outputDir.resolve("pivot");
    systemOut = tapSystemOut(() -> CalcAlleleFrequencies.main(new String[] {
        "-i", outputDir.toString(),
        "-o", pivotOutputDir.toString(),
        "-pc", "18"
    }));
    System.out.println(systemOut);
  }
}

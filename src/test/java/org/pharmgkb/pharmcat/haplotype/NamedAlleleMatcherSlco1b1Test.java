package org.pharmgkb.pharmcat.haplotype;

import java.nio.file.Path;
import java.util.List;
import com.google.common.collect.Lists;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.TestInfo;
import org.pharmgkb.common.util.PathUtils;
import org.pharmgkb.pharmcat.TestUtils;
import org.pharmgkb.pharmcat.haplotype.model.Result;
import org.pharmgkb.pharmcat.util.DataManager;

import static org.pharmgkb.pharmcat.haplotype.NamedAlleleMatcherTest.assertDiplotypePairs;
import static org.pharmgkb.pharmcat.haplotype.NamedAlleleMatcherTest.testMatchNamedAlleles;


/**
 *  * Test calling for SLCO1B1.
 *
 * @author Lester Carter
 */
public class NamedAlleleMatcherSlco1b1Test {
  private static final Path sf_definitionFile = DataManager.DEFAULT_DEFINITION_DIR.resolve("SLCO1B1_translation.json");


  @BeforeAll
  static void prepare() {
    //TestUtils.setSaveTestOutput(true);
  }

  @AfterEach
  void deleteDirectory(TestInfo testInfo) {
    TestUtils.deleteTestOutputDirectory(testInfo);
  }


  @Test
  void slco1b1s1s1() throws Exception {
    // Test *1/*1

    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/SLCO1B1/s1s1.vcf");
    List<String> expectedMatches = Lists.newArrayList("*1/*1");

    Result result = testMatchNamedAlleles(sf_definitionFile, vcfFile);
    assertDiplotypePairs(expectedMatches, result);
  }

  @Test
  void slco1b1s5s15() throws Exception {
    // Test *5/*15

    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/SLCO1B1/s5s15.vcf");
    List<String> expectedMatches = Lists.newArrayList("*5/*15");

    Result result = testMatchNamedAlleles(sf_definitionFile, vcfFile);
    assertDiplotypePairs(expectedMatches, result);
  }

  @Test
  void slco1b1s1s15() throws Exception {
    // Test *1/*15. Except we can't distinguish *1B/*5.

    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/SLCO1B1/s1s15.vcf");
    List<String> expectedMatches = Lists.newArrayList("*1/*15","*5/*37");

    Result result = testMatchNamedAlleles(sf_definitionFile, vcfFile);
    assertDiplotypePairs(expectedMatches, result);
  }
}

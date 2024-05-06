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
 * Test calling for CYP3A5.
 *
 * @author Lester Carter
 */
public class NamedAlleleMatcherCyp3a5Test {
  private static final Path sf_definitionFile = DataManager.DEFAULT_DEFINITION_DIR.resolve("CYP3A5_translation.json");


  @BeforeAll
  static void prepare() {
    //TestUtils.setSaveTestOutput(true);
  }

  @AfterEach
  void deleteDirectory(TestInfo testInfo) {
    TestUtils.deleteTestOutputDirectory(testInfo);
  }


  @Test
  void cyp3a5s3s9() throws Exception {
    // Test *3/*9.

    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/cyp3a5/s3s9.vcf");
    List<String> expectedMatches = Lists.newArrayList("*3/*9");

    Result result = testMatchNamedAlleles(sf_definitionFile, vcfFile);
    assertDiplotypePairs(expectedMatches, result);
  }

  @Test
  void cyp3a5s1s7() throws Exception {
    // Test *1/*7.  Het in rs41303343 position, a single insertion.

    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/cyp3a5/s1s7.vcf");
    List<String> expectedMatches = Lists.newArrayList("*1/*7");

    Result result = testMatchNamedAlleles(sf_definitionFile, vcfFile);
    assertDiplotypePairs(expectedMatches, result);
  }


  @Test
  void cyp3a5s3s9Homozygous() throws Exception {
    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/cyp3a5/s3s9-homozygous.vcf");
    List<String> expectedMatches = Lists.newArrayList("*3/*9");

    Result result = testMatchNamedAlleles(sf_definitionFile, vcfFile);
    assertDiplotypePairs(expectedMatches, result);
  }


  /**
   * Test example analogous to second scoring example in
   * <a href="https://pharmcat.org/methods/NamedAlleleMatcher-101">https://pharmcat.org/methods/NamedAlleleMatcher-101/</a>
   * but with *5 (in this case *7 for cyp3a5 missing).
   */
  @Test
  void s1s7missing() throws Exception {
    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/cyp3a5/s1s7missing.vcf");

    // just top match
    Result result = testMatchNamedAlleles(sf_definitionFile, vcfFile);
    List<String> expectedMatches = Lists.newArrayList("*1/*1");
    assertDiplotypePairs(expectedMatches, result);

    // All matches
    // TODO(markwoon): this one does NOT assume default
    Result result2 = testMatchNamedAlleles(sf_definitionFile, vcfFile, false, false, false);
    List<String> expectedMatches2 = Lists.newArrayList("*1/*1");
    assertDiplotypePairs(expectedMatches2, result2);

    // Don't presume reference and take all hits
    Result result3 = testMatchNamedAlleles(sf_definitionFile, vcfFile, false, false, false);
    List<String> expectedMatches3 = Lists.newArrayList("*1/*1");
    assertDiplotypePairs(expectedMatches3, result3);

  }

  @Test
  void s1s1rs776746missing() throws Exception {
    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/cyp3a5/s1s1rs776746missing.vcf");

    // TODO(markwoon): this one does NOT assume default
    Result result2 = testMatchNamedAlleles(sf_definitionFile, vcfFile, true, false, false);
    List<String> expectedMatches2 = Lists.newArrayList("*1/*1");
    assertDiplotypePairs(expectedMatches2, result2);
  }




}

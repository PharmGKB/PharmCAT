package org.pharmgkb.pharmcat.haplotype;

import java.nio.file.Path;
import java.util.List;
import com.google.common.collect.Lists;
import org.junit.Before;
import org.junit.Test;
import org.pharmgkb.common.util.PathUtils;
import org.pharmgkb.pharmcat.haplotype.model.Result;
import org.pharmgkb.pharmcat.util.DataManager;

import static org.pharmgkb.pharmcat.haplotype.NamedAlleleMatcherTest.assertDiplotypePairs;
import static org.pharmgkb.pharmcat.haplotype.NamedAlleleMatcherTest.testMatchNamedAlleles;


/**
 * JUnit test for {@link NamedAlleleMatcher#callDiplotypes(MatchData, boolean)}.
 *
 * @author Lester Carter
 */
public class NamedAlleleMatcherCyp3a5Test {
  private Path m_definitionFile;

  @Before
  public void before() throws Exception {
    m_definitionFile = DataManager.DEFAULT_DEFINITION_DIR.resolve("CYP3A5_translation.json");
  }

  @Test
  public void cyp3a5s3s9() throws Exception {
    // Test *3/*9.

    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/cyp3a5/s3s9.vcf");
    List<String> expectedMatches = Lists.newArrayList("*3/*9");

    Result result = testMatchNamedAlleles(m_definitionFile, vcfFile);
    assertDiplotypePairs(expectedMatches, result);
  }

  @Test
  public void cyp3a5s1s7() throws Exception {
    // Test *1/*7.  Het in rs41303343 position, a single insertion.

    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/cyp3a5/s1s7.vcf");
    List<String> expectedMatches = Lists.newArrayList("*1/*7");

    Result result = testMatchNamedAlleles(m_definitionFile, vcfFile);
    assertDiplotypePairs(expectedMatches, result);
  }


  @Test
  public void cyp3a5s3s9Homozygous() throws Exception {
    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/cyp3a5/s3s9-homozygous.vcf");
    List<String> expectedMatches = Lists.newArrayList("*1/*1");

    Result result = testMatchNamedAlleles(m_definitionFile, vcfFile);
    assertDiplotypePairs(expectedMatches, result);
  }

  //org/pharmgkb/pharmcat/haplotype/cyp3a5/s1s7missing.vcf

  @Test
  /* Test example analogous to second scoring example in
  https://github.com/PharmGKB/PharmCAT/wiki/NamedAlleleMatcher-101
  but with *5 (in this case *7 for cyp3a5 missing
  */

  public void s1s7missing() throws Exception {
    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/cyp3a5/s1s7missing.vcf");

    // just top match
    Result result = testMatchNamedAlleles(m_definitionFile, vcfFile);
    List<String> expectedMatches = Lists.newArrayList("*1/*1");
    assertDiplotypePairs(expectedMatches, result);

    // All matches
    Result result2 = testMatchNamedAlleles(m_definitionFile, vcfFile,false, false, false, false);
    List<String> expectedMatches2 = Lists.newArrayList("*1/*1");
    assertDiplotypePairs(expectedMatches2, result2);

    // Don't presume reference and take all hits
    Result result3 = testMatchNamedAlleles(m_definitionFile, vcfFile, true, false, false, false);
    List<String> expectedMatches3 = Lists.newArrayList("*1/*1");
    assertDiplotypePairs(expectedMatches3, result3);

  }

  @Test
  public void s1s1rs776746missing() throws Exception {
    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/cyp3a5/s1s1rs776746missing.vcf");

    Result result2 = testMatchNamedAlleles(m_definitionFile, vcfFile,false, true, false, false);
    List<String> expectedMatches2 = Lists.newArrayList("*1/*1");
    assertDiplotypePairs(expectedMatches2, result2);
  }




}

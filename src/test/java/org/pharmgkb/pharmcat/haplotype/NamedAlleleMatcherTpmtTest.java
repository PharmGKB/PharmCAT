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
public class NamedAlleleMatcherTpmtTest {
  private Path m_definitionFile;


  @Before
  public void before() throws Exception {
    m_definitionFile = DataManager.DEFAULT_DEFINITION_DIR.resolve("TPMT_translation.json");
  }


  @Test
  public void tpmts1s1() throws Exception {
    // Test *1/*1
    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/TPMT/s1s1.vcf");
    List<String> expectedMatches = Lists.newArrayList("*1/*1");

    Result result = testMatchNamedAlleles(m_definitionFile, vcfFile);
    assertDiplotypePairs(expectedMatches, result);
  }

  @Test
  public void tpmts1s1s() throws Exception {
    // without exemptions
    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/TPMT/s1s1s.vcf");
    List<String> expectedMatches = Lists.newArrayList("*1S/*1S");

    Result result = NamedAlleleMatcherTest.testMatchNamedAlleles(m_definitionFile, vcfFile, true, false, true, false);
    assertDiplotypePairs(expectedMatches, result);

    // with exemptions
    expectedMatches = Lists.newArrayList();

    result = NamedAlleleMatcherTest.testMatchNamedAlleles(m_definitionFile, vcfFile, true, false, true, true);
    assertDiplotypePairs(expectedMatches, result);
  }

  @Test
  public void tpmts3bs3c() throws Exception {
    // Test *3b/*3c.  However due to lack of phasing *1/*3a is also an option
    // However we are only taking top option at the moment, so *1/*3A (which scores higher) is valid single return
    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/TPMT/s3bs3c.vcf");
    List<String> expectedMatches = Lists.newArrayList("*1/*3A");

    Result result = testMatchNamedAlleles(m_definitionFile, vcfFile);
    assertDiplotypePairs(expectedMatches, result);
  }

  @Test
  public void tpmts3as3b() throws Exception {
    // Test *3a/*3b.
    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/TPMT/s3as3b.vcf");
    List<String> expectedMatches = Lists.newArrayList("*3A/*3B");

    Result result = testMatchNamedAlleles(m_definitionFile, vcfFile);
    assertDiplotypePairs(expectedMatches, result);
  }

  @Test
  public void tpmts16s22() throws Exception {
    // Test *16/*22. rs144041067 is bi-allelic
    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/TPMT/s16s22.vcf");
    List<String> expectedMatches = Lists.newArrayList("*16/*22");

    Result result = testMatchNamedAlleles(m_definitionFile, vcfFile);
    assertDiplotypePairs(expectedMatches, result);
  }

  /**
   * Test to make sure homo *1S alleles are ignored and that other present alleles will be called properly
   */
  @Test
  public void tpmtStar1SHomo3AHet() throws Exception {
    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/TPMT/s1ss1ss3.vcf");
    List<String> expectedMatches = Lists.newArrayList("*1/*3A");

    Result result = testMatchNamedAlleles(m_definitionFile, vcfFile);
    assertDiplotypePairs(expectedMatches, result);
  }
}

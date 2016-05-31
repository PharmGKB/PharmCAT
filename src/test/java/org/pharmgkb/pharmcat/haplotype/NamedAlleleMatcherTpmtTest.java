package org.pharmgkb.pharmcat.haplotype;

import java.nio.file.Path;
import java.util.List;
import com.google.common.collect.Lists;
import org.junit.Before;
import org.junit.Test;
import org.pharmgkb.common.util.PathUtils;
import org.pharmgkb.pharmcat.haplotype.model.Result;

import static org.pharmgkb.pharmcat.haplotype.NamedAlleleMatcherTest.assertDiplotypePairs;
import static org.pharmgkb.pharmcat.haplotype.NamedAlleleMatcherTest.testMatchNamedAlleles;


/**
 * JUnit test for {@link NamedAlleleMatcher#callDiplotypes(MatchData)}.
 *
 * @author Lester Carter
 */
public class NamedAlleleMatcherTpmtTest {
  private Path m_definitionFile;


  @Before
  public void before() throws Exception {
    m_definitionFile =  PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/TPMT_translation.json");
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
  public void tpmts3bs3c() throws Exception {
    // Test *3b/*3c.  However due to lack of phasing *1/*3a is also an option
    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/TPMT/s3bs3c.vcf");
    List<String> expectedMatches = Lists.newArrayList("*1/*3A","*3B/*3C");

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
}

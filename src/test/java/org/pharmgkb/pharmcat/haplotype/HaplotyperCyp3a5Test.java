package org.pharmgkb.pharmcat.haplotype;

import java.nio.file.Path;
import java.util.List;
import com.google.common.collect.Lists;
import org.junit.Before;
import org.junit.Test;
import org.pharmgkb.pharmcat.TestUtil;
import org.pharmgkb.pharmcat.haplotype.model.HaplotyperResult;

import static org.pharmgkb.pharmcat.haplotype.HaplotyperTest.assertDiplotypePairs;
import static org.pharmgkb.pharmcat.haplotype.HaplotyperTest.testCallHaplotype;


/**
 * JUnit test for {@link Haplotyper#callDiplotypes(MatchData)}.
 *
 * @author Lester Carter
 */
public class HaplotyperCyp3a5Test {
  private Path m_definitionFile;

  @Before
  public void before() throws Exception {
    m_definitionFile =  TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/CYP3A5_translation.json");
  }

  @Test
  public void cyp3a5s3d9() throws Exception {
    // Test *3/*9.  Note Y in 99672916 position

    Path vcfFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/cyp3a5/s3s9.vcf");
    List<String> expectedMatches = Lists.newArrayList("*3/*9");

    HaplotyperResult result = testCallHaplotype(m_definitionFile, vcfFile);
    assertDiplotypePairs(expectedMatches, result);
  }

  @Test
  public void cyp3a5s1s7() throws Exception {
    // Test *1/*7.  Het in rs41303343 position, a single insertion.

    Path vcfFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/cyp3a5/s1s7.vcf");
    List<String> expectedMatches = Lists.newArrayList("*1/*7");

    HaplotyperResult result = testCallHaplotype(m_definitionFile, vcfFile);
    assertDiplotypePairs(expectedMatches, result);
  }
}

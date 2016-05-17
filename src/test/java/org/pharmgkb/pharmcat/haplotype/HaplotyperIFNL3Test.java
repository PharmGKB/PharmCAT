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
public class HaplotyperIFNL3Test {
  private Path m_definitionFile;

  @Before
  public void before() throws Exception {
    m_definitionFile =  TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/IFNL3_translation.json");
  }

  @Test
  public void rs12979860CC() throws Exception {
    // Test rs12979860 CC

    Path vcfFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/IFNL3/rs12979860CC.vcf");
    List<String> expectedMatches = Lists.newArrayList("rs12979860C/rs12979860C");

    HaplotyperResult result = testCallHaplotype(m_definitionFile, vcfFile);
    assertDiplotypePairs(expectedMatches, result);
  }

  @Test
  public void rs12979860CT() throws Exception {
    // Test rs12979860 CC

    Path vcfFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/IFNL3/rs12979860CT.vcf");
    List<String> expectedMatches = Lists.newArrayList("rs12979860C/rs12979860T");

    HaplotyperResult result = testCallHaplotype(m_definitionFile, vcfFile);
    assertDiplotypePairs(expectedMatches, result);
  }

  @Test
  public void rs12979860TT() throws Exception {
    // Test rs12979860 CC

    Path vcfFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/IFNL3/rs12979860TT.vcf");
    List<String> expectedMatches = Lists.newArrayList("rs12979860T/rs12979860T");

    HaplotyperResult result = testCallHaplotype(m_definitionFile, vcfFile);
    assertDiplotypePairs(expectedMatches, result);
  }
}

package org.pharmgkb.pharmcat.haplotype;

import java.nio.file.Path;
import java.util.List;
import com.google.common.collect.Lists;
import org.junit.Before;
import org.junit.Test;
import org.pharmgkb.pharmcat.TestUtil;
import org.pharmgkb.pharmcat.haplotype.model.DiplotypeMatch;


/**
 * JUnit test for {@link Haplotyper#callDiplotypes(MatchData)}.
 *
 * @author Lester Carter
 */
public class HaplotyperVkorc1Test {
  private Path m_jsonFile;

  @Before
  public void before() throws Exception {
    m_jsonFile =  TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/VKORC1_translation.json");
  }

  @Test
  public void vkorc1gg() throws Exception {
    // Test -1639G/-1639G

    Path vcfFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/VKORC1/-1639G-1639G.vcf");
    List<DiplotypeMatch> matches = HaplotyperTest.testCallHaplotype(m_jsonFile, vcfFile);

    List<String> expectedMatches = Lists.newArrayList("-1639G/-1639G");
    TestUtil.assertDiplotypePairs(expectedMatches, matches);
  }

  @Test
  public void vkorc1ga() throws Exception {
    // Test -1639G/-1639A

    Path vcfFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/VKORC1/-1639G-1639A.vcf");
    List<DiplotypeMatch> matches = HaplotyperTest.testCallHaplotype(m_jsonFile, vcfFile);

    List<String> expectedMatches = Lists.newArrayList("-1639A/-1639G");
    TestUtil.assertDiplotypePairs(expectedMatches, matches);
  }

  @Test
  public void vkorc1aa() throws Exception {
    // Test -1639A/-1639A

    Path vcfFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/VKORC1/-1639A-1639A.vcf");
    List<DiplotypeMatch> matches = HaplotyperTest.testCallHaplotype(m_jsonFile, vcfFile);

    List<String> expectedMatches = Lists.newArrayList("-1639A/-1639A");
    TestUtil.assertDiplotypePairs(expectedMatches, matches);
  }


}

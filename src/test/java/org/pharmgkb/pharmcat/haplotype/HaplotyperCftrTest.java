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
public class HaplotyperCftrTest {
  private Path m_jsonFile;

  @Before
  public void before() throws Exception {
    m_jsonFile =  TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/CFTR_translation.json");
  }

  // TODO: Lester - need to double check vcf.

  @Test
  public void cftrF508delF508del() throws Exception {
    // Test *1/*1 TODO: Lester - check that the star is defined correctly. Fails now as *1 contains ATATATATATATATAA
    Path vcfFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/CFTR/F508delF508del.vcf");
    List<DiplotypeMatch> matches = HaplotyperTest.testCallHaplotype(m_jsonFile, vcfFile);

    List<String> expectedMatches = Lists.newArrayList("F508del/F508del");
    TestUtil.assertDiplotypePairs(expectedMatches, matches);
  }

}

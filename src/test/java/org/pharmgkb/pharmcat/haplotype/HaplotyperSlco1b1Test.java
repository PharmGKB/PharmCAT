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
public class HaplotyperSlco1b1Test {
  private Path m_tsvFile;

  @Before
  public void before() throws Exception {
    m_tsvFile =  TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/SLCO1B1_translation.json");
  }

  @Test
  public void slco1b1s1as1a() throws Exception {
    // Test *1/*1

    Path vcfFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/SLCO1B1/s1as1a.vcf");
    List<DiplotypeMatch> matches = HaplotyperTest.testCallHaplotype(m_tsvFile, vcfFile);

    List<String> expectedMatches = Lists.newArrayList("*1A/*1A");
    TestUtil.assertDiplotypePairs(expectedMatches, matches);
  }

  @Test
  public void slco1b1s5s15() throws Exception {
    // Test *5/*15

    Path vcfFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/SLCO1B1/s5s15.vcf");
    List<DiplotypeMatch> matches = HaplotyperTest.testCallHaplotype(m_tsvFile, vcfFile);

    List<String> expectedMatches = Lists.newArrayList("*5/*15");
    TestUtil.assertDiplotypePairs(expectedMatches, matches);
  }

  @Test
  public void slco1b1s1bs15() throws Exception {
    // Test *1b/*15. Except we can't distinguish *1B/*5 or *5/*29

    Path vcfFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/SLCO1B1/s1bs15.vcf");
    List<DiplotypeMatch> matches = HaplotyperTest.testCallHaplotype(m_tsvFile, vcfFile);

    List<String> expectedMatches = Lists.newArrayList("*1A/*15","*1B/*5","*5/*29");
    TestUtil.assertDiplotypePairs(expectedMatches, matches);
  }

  @Test
  public void slco1b17s1bs21() throws Exception {
    // Test *17/*21

    Path vcfFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/SLCO1B1/s17s21.vcf");
    List<DiplotypeMatch> matches = HaplotyperTest.testCallHaplotype(m_tsvFile, vcfFile);

    List<String> expectedMatches = Lists.newArrayList("*17/*21");
    TestUtil.assertDiplotypePairs(expectedMatches, matches);
  }



}

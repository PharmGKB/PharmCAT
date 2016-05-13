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
public class HaplotyperCyp2c9Test {
  private Path m_tsvFile;

  @Before
  public void before() throws Exception {
    m_tsvFile =  TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/CYP2C9_translation.json");
  }

  @Test
  public void cyp32c9s1s1() throws Exception {
    // Test *1/*1

    Path vcfFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/cyp2c9/s1s1.vcf");
    List<DiplotypeMatch> matches = HaplotyperTest.testCallHaplotype(m_tsvFile, vcfFile);

    List<String> expectedMatches = Lists.newArrayList("*1/*1");
    TestUtil.assertDiplotypePairs(expectedMatches, matches);
  }

  @Test
  public void cyp32c9s2s3() throws Exception {
    // Test *2/*3

    Path vcfFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/cyp2c9/s2s3.vcf");
    List<DiplotypeMatch> matches = HaplotyperTest.testCallHaplotype(m_tsvFile, vcfFile);

    List<String> expectedMatches = Lists.newArrayList("*2/*3");
    TestUtil.assertDiplotypePairs(expectedMatches, matches);
  }

  @Test
  public void cyp32c9s2s24() throws Exception {
    // Test *2/*24, but also matches *1/*24

    /*
    TODO - fails.  But says it doesn't!:

    *1/*24 (108), *2/*24 (108)
    Printing to /Users/lester/IdeaProjects/cpic-annotator/build/resources/test/org/pharmgkb/pharmcat/haplotype/cyp2c9/s2s24.html
    Expected: [*1/*24, *2/*24]
    Got:      [*1/*24, *2/*24]

    java.lang.AssertionError: Did not get expected matches

     */


    Path vcfFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/cyp2c9/s2s24.vcf");
    List<DiplotypeMatch> matches = HaplotyperTest.testCallHaplotype(m_tsvFile, vcfFile);

    List<String> expectedMatches = Lists.newArrayList("*1/*24", "*2/*24");
    TestUtil.assertDiplotypePairs(expectedMatches, matches);
  }

  @Test
  public void cyp32c9s2s24hom() throws Exception {
    // Test *2/*24 with hom at rs749060448

    /* TODO - Lester: Fails. Micelle expected *2/*24. We are checking logic on our side now.

    *24/*24 (108)
    Printing to /Users/lester/IdeaProjects/cpic-annotator/build/resources/test/org/pharmgkb/pharmcat/haplotype/cyp2c9/s2s24hom.html
    Expected: [*2/*24]
    Got:      [*24/*24]
     */

    Path vcfFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/cyp2c9/s2s24hom.vcf");
    List<DiplotypeMatch> matches = HaplotyperTest.testCallHaplotype(m_tsvFile, vcfFile, true, false, true);

    List<String> expectedMatches = Lists.newArrayList("*2/*24");
    TestUtil.assertDiplotypePairs(expectedMatches, matches);
  }
}

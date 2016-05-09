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
public class HaplotyperUgt1a1Test {
  private Path m_jsonFile;

  // TODO: lester - work in progress

  @Before
  public void before() throws Exception {
     m_jsonFile =  TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/UGT1A1_translation.json");
  }

  @Test
  public void ugt1a1s1s1() throws Exception {
    // Test *1/*1 TODO: Lester - check that the star is defined correctly *1 contains ATATATATATATATAA
    Path vcfFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/UGT1A1/s1s1.vcf");
    List<DiplotypeMatch> matches = HaplotyperTest.testCallHaplotype( m_jsonFile, vcfFile);

    List<String> expectedMatches = Lists.newArrayList("*1/*1");
    TestUtil.assertDiplotypePairs(expectedMatches, matches);
  }

  @Test
  public void ugt1a1s28s37() throws Exception {
    // Test *28/*37 contains TA repeat
    Path vcfFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/UGT1A1/s28s37.vcf");
    List<DiplotypeMatch> matches = HaplotyperTest.testCallHaplotype( m_jsonFile, vcfFile);

    List<String> expectedMatches = Lists.newArrayList("*28/*37");
    TestUtil.assertDiplotypePairs(expectedMatches, matches);
  }

  @Test
  public void ugt1a1s28s80() throws Exception {
    // Test *28/*80 Another TA repeat example
    Path vcfFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/UGT1A1/s28s80.vcf");
    List<DiplotypeMatch> matches = HaplotyperTest.testCallHaplotype( m_jsonFile, vcfFile);

    List<String> expectedMatches = Lists.newArrayList("*28/*80");
    TestUtil.assertDiplotypePairs(expectedMatches, matches);
  }

  @Test
  public void ugt1a1s6s6() throws Exception {
    // Test *6/*6
    Path vcfFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/UGT1A1/s6s6.vcf");
    List<DiplotypeMatch> matches = HaplotyperTest.testCallHaplotype( m_jsonFile, vcfFile);

    List<String> expectedMatches = Lists.newArrayList("*6/*6");
    TestUtil.assertDiplotypePairs(expectedMatches, matches);
  }

}

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
 * @author Mark Woon
 */
public class HaplotyperCyp2c19Test {
  private Path m_definitionFile;

  @Before
  public void before() throws Exception {
     m_definitionFile =  TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/CYP2C19_translation.json");
  }


  @Test
  public void cyp2c19s1s1() throws Exception {
    // Simple *1/*1 reference test
    Path vcfFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/cyp2c19/s1s1.vcf");
    List<String> expectedMatches = Lists.newArrayList("*1/*1");

    HaplotyperResult result = testCallHaplotype(m_definitionFile, vcfFile);
    assertDiplotypePairs(expectedMatches, result);
  }


  @Test
  public void cyp2c19s1s2() throws Exception {
    // Test simple case of one heterozygous snp
    Path vcfFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/cyp2c19/s1s2.vcf");
    List<String> expectedMatches = Lists.newArrayList("*1/*2");

    HaplotyperResult result = testCallHaplotype(m_definitionFile, vcfFile);
    assertDiplotypePairs(expectedMatches, result);
  }


  @Test
  public void cyp2c19s1s4b() throws Exception {
    // Test *1 *4b. s1s4b-longest wins(not s4as17)
    Path vcfFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/cyp2c19/s1s4b.vcf");
    List<String> expectedMatches = Lists.newArrayList("*1/*4B", "*4A/*17");

    HaplotyperResult result = testCallHaplotype(m_definitionFile, vcfFile);
    assertDiplotypePairs(expectedMatches, result);
  }


  @Test
  public void cyp2c19s1s4bMissingCalls() throws Exception {
    // Test *1 *4b - but with some of the calls deleted from the vcf.  As requested by Michelle.
    Path vcfFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/cyp2c19/s1s4bMissing.vcf");
    List<String> expectedMatches = Lists.newArrayList("*1/*4B", "*4A/*17");

    HaplotyperResult result = testCallHaplotype(m_definitionFile, vcfFile);
    assertDiplotypePairs(expectedMatches, result);
  }


  @Test
  public void cyp2c19s1s28() throws Exception {
    // Test *1 *28. The longest possible match should win.
    Path vcfFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/cyp2c19/s1s28.vcf");
    List<String> expectedMatches = Lists.newArrayList("*1/*28");

    HaplotyperResult result = testCallHaplotype(m_definitionFile, vcfFile);
    assertDiplotypePairs(expectedMatches, result);
  }


  @Test
  public void cyp2c19s2s2() throws Exception {
    // Test simple case of one homozygous snp
    Path vcfFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/cyp2c19/s2s2.vcf");
    List<String> expectedMatches = Lists.newArrayList("*2/*2");

    HaplotyperResult result = testCallHaplotype(m_definitionFile, vcfFile);
    assertDiplotypePairs(expectedMatches, result);
  }


  @Test
  public void cyp2c19s2s3() throws Exception {
    // Test case of two snps - *2 and *3
    Path vcfFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/cyp2c19/s2s3.vcf");
    List<String> expectedMatches = Lists.newArrayList("*2/*3");

    HaplotyperResult result = testCallHaplotype(m_definitionFile, vcfFile);
    assertDiplotypePairs(expectedMatches, result);
  }


  @Test
  public void cyp2c19s4as4b() throws Exception {
    // Test *4a *4b. het and homo rsids
    Path vcfFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/cyp2c19/s4as4b.vcf");
    List<String> expectedMatches = Lists.newArrayList("*4A/*4B");

    HaplotyperResult result = testCallHaplotype(m_definitionFile, vcfFile);
    assertDiplotypePairs(expectedMatches, result);
  }


  @Test
  public void cyp2c19s4bs17() throws Exception {
    // Test *4b/*17
    Path vcfFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/cyp2c19/s4bs17.vcf");
    List<String> expectedMatches = Lists.newArrayList("*4B/*17");

    HaplotyperResult result = testCallHaplotype(m_definitionFile, vcfFile);
    assertDiplotypePairs(expectedMatches, result);
  }


  @Test
  public void cyp2c19s15s28() throws Exception {
    // Test *15 *28. The shared position is homo
    Path vcfFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/cyp2c19/s15s28.vcf");
    List<String> expectedMatches = Lists.newArrayList("*15/*28");

    HaplotyperResult result = testCallHaplotype(m_definitionFile, vcfFile);
    assertDiplotypePairs(expectedMatches, result);
  }


  @Test
  public void cyp2c19sUnks17() throws Exception {
    // Test *Unk/*17 - only one haplotype matches, so no diploid match
    Path vcfFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/cyp2c19/sUnks17.vcf");
    // no matches
    List<String> expectedMatches = Lists.newArrayList();

    HaplotyperResult result = testCallHaplotype(m_definitionFile, vcfFile);
    assertDiplotypePairs(expectedMatches, result);
  }
}

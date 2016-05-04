package org.pharmgkb.pharmcat.haplotype;

import java.nio.file.Path;
import java.util.List;
import java.util.SortedMap;
import com.google.common.collect.Lists;
import org.junit.Before;
import org.junit.Test;
import org.pharmgkb.pharmcat.TestUtil;
import org.pharmgkb.pharmcat.haplotype.model.DiplotypeMatch;


/**
 * JUnit test for {@link Haplotyper#callDiplotypes(SortedMap, String)}.
 *
 * @author Lester Carter
 */
public class HaplotyperTpmtTest {
  private Path m_tsvFile;

  @Before
  public void before() throws Exception {
    m_tsvFile =  TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/TPMT.tsv");
  }


  @Test
  public void tpmts1s1() throws Exception {
    // Test *1/*1
    Path vcfFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/TPMT/s1s1.vcf");
    List<DiplotypeMatch> matches = HaplotyperTest.testCallHaplotype(m_tsvFile, vcfFile);

    List<String> expectedMatches = Lists.newArrayList("*1/*1");
    TestUtil.assertDiplotypePairs(expectedMatches, matches);
  }

  @Test
  public void tpmts3bs3c() throws Exception {
    // Test *3b/*3c.  However due to lack of phasing *1/*3a is also an option
    Path vcfFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/TPMT/s3bs3c.vcf");
    List<DiplotypeMatch> matches = HaplotyperTest.testCallHaplotype(m_tsvFile, vcfFile);

    List<String> expectedMatches = Lists.newArrayList("*1/*3A","*3B/*3C");
    TestUtil.assertDiplotypePairs(expectedMatches, matches);
  }

  @Test
  public void tpmts3as3b() throws Exception {
    // Test *3a/*3b.
    Path vcfFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/TPMT/s3as3b.vcf");
    List<DiplotypeMatch> matches = HaplotyperTest.testCallHaplotype(m_tsvFile, vcfFile);

    List<String> expectedMatches = Lists.newArrayList("*3A/*3B");
    TestUtil.assertDiplotypePairs(expectedMatches, matches);
  }

  @Test
  public void tpmts16s22() throws Exception {
    // Test *16/*22. rs144041067 is bi-allelic
    Path vcfFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/TPMT/s16s22.vcf");
    List<DiplotypeMatch> matches = HaplotyperTest.testCallHaplotype(m_tsvFile, vcfFile);

    List<String> expectedMatches = Lists.newArrayList("*16/*22");
    TestUtil.assertDiplotypePairs(expectedMatches, matches);
  }



}

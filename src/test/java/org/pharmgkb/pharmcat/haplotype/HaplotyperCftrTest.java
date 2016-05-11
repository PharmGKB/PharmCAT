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


  @Test
  public void cftrReferenceReference() throws Exception {
    // Test reference
    Path vcfFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/CFTR/refref.vcf");
    List<DiplotypeMatch> matches = HaplotyperTest.testCallHaplotype(m_jsonFile, vcfFile);

    List<String> expectedMatches = Lists.newArrayList("Reference/Reference");
    TestUtil.assertDiplotypePairs(expectedMatches, matches);
  }


  @Test
  public void cftrF508delF508del() throws Exception {
    // Test F508del(TCT)/F508del(TCT)
    Path vcfFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/CFTR/F508delF508del.vcf");
    List<DiplotypeMatch> matches = HaplotyperTest.testCallHaplotype(m_jsonFile, vcfFile);

    List<String> expectedMatches = Lists.newArrayList("F508del(TCT)/F508del(TCT)");
    TestUtil.assertDiplotypePairs(expectedMatches, matches);
  }

  @Test
  public void cftrWt5T() throws Exception {
    // Test Reference/5T - TODO - assumption is that repeat is represented correctly
    Path vcfFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/CFTR/ref5t.vcf");
    List<DiplotypeMatch> matches = HaplotyperTest.testCallHaplotype(m_jsonFile, vcfFile);

    List<String> expectedMatches = Lists.newArrayList("5T/Reference");
    TestUtil.assertDiplotypePairs(expectedMatches, matches);
  }
}

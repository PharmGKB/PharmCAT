package org.pharmgkb.pharmcat.haplotype;

import java.nio.file.Path;
import java.util.List;
import com.google.common.collect.Lists;
import org.junit.Before;
import org.junit.Test;
import org.pharmgkb.common.util.PathUtils;
import org.pharmgkb.pharmcat.haplotype.model.HaplotyperResult;

import static org.pharmgkb.pharmcat.haplotype.HaplotyperTest.assertDiplotypePairs;
import static org.pharmgkb.pharmcat.haplotype.HaplotyperTest.testCallHaplotype;


/**
 * JUnit test for {@link Haplotyper#callDiplotypes(MatchData)}.
 *
 * @author Lester Carter
 */
public class HaplotyperUgt1a1Test {
  private Path m_definitionFile;

  // TODO: lester - work in progress.  This needs checking against real world examples.

  @Before
  public void before() throws Exception {
    m_definitionFile =  PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/UGT1A1_translation.json");
  }

  @Test
  public void ugt1a1s1s1() throws Exception {
    // Test *1/*1
    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/UGT1A1/s1s1.vcf");
    List<String> expectedMatches = Lists.newArrayList("*1/*1");

    HaplotyperResult result = testCallHaplotype(m_definitionFile, vcfFile);
    assertDiplotypePairs(expectedMatches, result);
  }

  @Test
  public void ugt1a1s28s37() throws Exception {
    // Test *28/*37 contains TA repeat
    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/UGT1A1/s28s37.vcf");
    List<String> expectedMatches = Lists.newArrayList("*28/*37");

    HaplotyperResult result = testCallHaplotype(m_definitionFile, vcfFile);
    assertDiplotypePairs(expectedMatches, result);
  }

  //@Test
  public void ugt1a1s28s80() throws Exception {
    // Test *28/*80 Another TA repeat example  - as above. Skipping for now as same logic.
    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/UGT1A1/s28s80.vcf");
    List<String> expectedMatches = Lists.newArrayList("*28/*80");

    HaplotyperResult result = testCallHaplotype(m_definitionFile, vcfFile);
    assertDiplotypePairs(expectedMatches, result);
  }

  @Test
  public void ugt1a1s6s6() throws Exception {
    // Test *6/*6
    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/UGT1A1/s6s6.vcf");
    List<String> expectedMatches = Lists.newArrayList("*6/*6");

    HaplotyperResult result = testCallHaplotype(m_definitionFile, vcfFile);
    assertDiplotypePairs(expectedMatches, result);
  }
}

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
public class HaplotyperSlco1b1Test {
  private Path m_definitionFile;

  @Before
  public void before() throws Exception {
    m_definitionFile =  PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/SLCO1B1_translation.json");
  }

  @Test
  public void slco1b1s1as1a() throws Exception {
    // Test *1/*1

    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/SLCO1B1/s1as1a.vcf");
    List<String> expectedMatches = Lists.newArrayList("*1A/*1A");

    HaplotyperResult result = testCallHaplotype(m_definitionFile, vcfFile);
    assertDiplotypePairs(expectedMatches, result);
  }

  @Test
  public void slco1b1s5s15() throws Exception {
    // Test *5/*15

    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/SLCO1B1/s5s15.vcf");
    List<String> expectedMatches = Lists.newArrayList("*5/*15");

    HaplotyperResult result = testCallHaplotype(m_definitionFile, vcfFile);
    assertDiplotypePairs(expectedMatches, result);
  }

  @Test
  public void slco1b1s1as15() throws Exception {
    // Test *1a/*15. Except we can't distinguish *1B/*5.

    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/SLCO1B1/s1as15.vcf");
    List<String> expectedMatches = Lists.newArrayList("*1A/*15","*1B/*5");

    HaplotyperResult result = testCallHaplotype(m_definitionFile, vcfFile);
    assertDiplotypePairs(expectedMatches, result);
  }

  @Test
  public void slco1b1s1as15s1bs5Missing() throws Exception {
    /* Test *1a/*15. Except we can't distinguish *1B/*5.

    However in this case we are missing the final position in the file:
    chr12	21239158	rs140790673	C	T	.	PASS	assume-default	GT	0/0

    which forms part of the *29 definition, so we also get:
    *5/*29

    Output should report that *29 is only partially called.
     */

    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/SLCO1B1/s1as15s1bs5missing.vcf");
    List<String> expectedMatches = Lists.newArrayList("*1A/*15","*1B/*5", "*5/*29");

    HaplotyperResult result = testCallHaplotype(m_definitionFile, vcfFile);
    assertDiplotypePairs(expectedMatches, result);
  }

  @Test
  public void slco1b1s1as15s1bs5TwoMissing() throws Exception {
    /* Test two positions missing:
      chr12	21239145	rs200995543	C	T	.	PASS	only--star-3-4	GT	0/0
      chr12	21239158	rs140790673	C	T	.	PASS	second-star-29	GT	0/0


    *1a/*15. Except we can't distinguish *1B/*5. *5/*29 Matches as well,
     because the first position matches.  However this shows difference from
     single missing position. This time *34 should be reported as 'can't call'
     while *29 should be reported as ony partially called.

     TODO: test passes, but we need to make sure the output makes this distinction clear
     */

    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/SLCO1B1/s1as15s1bs5twomissing.vcf");
    List<String> expectedMatches = Lists.newArrayList("*1A/*15","*1B/*5", "*5/*29");

    HaplotyperResult result = testCallHaplotype(m_definitionFile, vcfFile);
    assertDiplotypePairs(expectedMatches, result);
  }


  @Test
  public void s1as15s1bs5NotCalled() throws Exception {
    // Test *1a/*15. Except we can't distinguish *1B/*5 or *5/*29, because final position is ./.

    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/SLCO1B1/s1as15s1bs5notcalled.vcf");
    // same as above, but using ./. to signify missing snp.
    List<String> expectedMatches = Lists.newArrayList("*1A/*15","*1B/*5", "*5/*29");

    HaplotyperResult result = testCallHaplotype(m_definitionFile, vcfFile);
    assertDiplotypePairs(expectedMatches, result);
  }


  @Test
  public void slco1b17s1bs21() throws Exception {
    // Test *17/*21

    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/SLCO1B1/s17s21.vcf");
    List<String> expectedMatches = Lists.newArrayList("*17/*21");

    HaplotyperResult result = testCallHaplotype(m_definitionFile, vcfFile);
    assertDiplotypePairs(expectedMatches, result);
  }
}

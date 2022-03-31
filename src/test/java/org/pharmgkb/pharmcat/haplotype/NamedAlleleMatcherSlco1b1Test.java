package org.pharmgkb.pharmcat.haplotype;

import java.nio.file.Path;
import java.util.List;
import com.google.common.collect.Lists;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.pharmgkb.common.util.PathUtils;
import org.pharmgkb.pharmcat.haplotype.model.Result;
import org.pharmgkb.pharmcat.util.DataManager;

import static org.pharmgkb.pharmcat.haplotype.NamedAlleleMatcherTest.assertDiplotypePairs;
import static org.pharmgkb.pharmcat.haplotype.NamedAlleleMatcherTest.testMatchNamedAlleles;


/**
 *  * Test calling for SLCO1B1.
 *
 * @author Lester Carter
 */
class NamedAlleleMatcherSlco1b1Test {
  private Path m_definitionFile;

  @BeforeEach
  void before() {
    m_definitionFile = DataManager.DEFAULT_DEFINITION_DIR.resolve("SLCO1B1_translation.json");
  }

  @Test
  void slco1b1s1s1() throws Exception {
    // Test *1/*1

    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/SLCO1B1/s1s1.vcf");
    List<String> expectedMatches = Lists.newArrayList("*1/*1");

    Result result = testMatchNamedAlleles(m_definitionFile, vcfFile);
    assertDiplotypePairs(expectedMatches, result);
  }

  @Test
  void slco1b1s5s15() throws Exception {
    // Test *5/*15

    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/SLCO1B1/s5s15.vcf");
    List<String> expectedMatches = Lists.newArrayList("*5/*15");

    Result result = testMatchNamedAlleles(m_definitionFile, vcfFile);
    assertDiplotypePairs(expectedMatches, result);
  }

  @Test
  void slco1b1s1s15() throws Exception {
    // Test *1/*15. Except we can't distinguish *1B/*5.

    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/SLCO1B1/s1s15.vcf");
    List<String> expectedMatches = Lists.newArrayList("*1/*15","*5/*37");

    Result result = testMatchNamedAlleles(m_definitionFile, vcfFile);
    assertDiplotypePairs(expectedMatches, result);
  }

  @Test
  void slco1b1s1as15s1bs5Missing() throws Exception {
    /* Test *1/*15. Except we can't distinguish *5/*37.

    However, in this case we are missing the final position in the file:
    chr12	21239158	rs140790673	C	T	.	PASS	assume-default	GT	0/0

    which forms part of the *29 definition, so we also get:
    *5/*29

    Output should report that *29 is only partially called.
     */

    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/SLCO1B1/s1as15s1bs5missing.vcf");
    List<String> expectedMatches = Lists.newArrayList("*1/*15", "*5/*29", "*5/*37");

    Result result = testMatchNamedAlleles(m_definitionFile, vcfFile);
    assertDiplotypePairs(expectedMatches, result);
  }

  @Test
  void slco1b1s1as15s1bs5TwoMissing() throws Exception {
    /* Test two positions missing:
      chr12	21239145	rs200995543	C	T	.	PASS	only--star-3-4	GT	0/0
      chr12	21239158	rs140790673	C	T	.	PASS	second-star-29	GT	0/0


    *1/*15. Except we can't distinguish *5/*37. *5/*29 Matches as well,
     because the first position matches.  However, this shows difference from
     single missing position. This time *34 should be reported as 'can't call'
     while *29 should be reported as ony partially called.
     */

    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/SLCO1B1/s1as15s1bs5twomissing.vcf");
    List<String> expectedMatches = Lists.newArrayList("*1/*15", "*5/*29", "*5/*37");

    Result result = testMatchNamedAlleles(m_definitionFile, vcfFile);
    assertDiplotypePairs(expectedMatches, result);
  }
}

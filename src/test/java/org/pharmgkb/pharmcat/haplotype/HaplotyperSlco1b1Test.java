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
  public void slco1b1s1as15() throws Exception {
    // Test *1a/*15. Except we can't distinguish *1B/*5.

    Path vcfFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/SLCO1B1/s1as15.vcf");
    List<DiplotypeMatch> matches = HaplotyperTest.testCallHaplotype(m_tsvFile, vcfFile);

    List<String> expectedMatches = Lists.newArrayList("*1A/*15","*1B/*5");
    TestUtil.assertDiplotypePairs(expectedMatches, matches);
  }

  @Test
  public void slco1b1s1as15s1bs5Missing() throws Exception {
    /* Test *1a/*15. Except we can't distinguish *1B/*5.

    However in this case we are missing the final position in the file:
    chr12	21239158	rs140790673	C	T	.	PASS	assume-default	GT	0/0

    which forms part of the *29 definition, so we also get:
    *5/*29

    At the moment the haplotyper removes from the tsv file analysis anything that is missing from the input file.
    In this case this changes the *29 definition to just be one snp, which then matches and hence we get the
    *5/*29 result as well.

    There a a few things that come out of this:

    1) If the normalized input vcf does not contain this position we currently consider this as a possibility.
    So we are showing the universe of all possible calls. If we do this we should have some way of showing
    in the report that this call is made using some unknown data.
    2) This means that the less calls we have the more potential diplotypes we can match.
    3) The input file vcf contains what we know we have.  Anything not in it we don't know.

    So therefore the normalized input vcf file should *always* contain all known positions, and this is the responsibility of the
    pre-parser. If it is not in the vcf it is presumed not read or not readable. Genotype "./." is also treated as missing.


    Basically this all comes down to a) deciding the missing data is a wild card (essentially what we do now by excluding it from the
    tsv analysis, b) the missing data is reference (which would remove the *29 possibility and make it behave more like I *think* Michelle
    expects) or c) Missing data really is missing and we can't call ANYTHING, even *1.

    It would be useful if the output html still showed all the tsv positions, as it's easier to scan along by eye.

    TLDR: we don't currently expect to match *5/*29. This needs discussing to make sure Michelle etc are okay adding this to the
    expected matches. The output does not show *29 was called using incomplete data.
     */

    Path vcfFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/SLCO1B1/s1as15s1bs5missing.vcf");
    List<DiplotypeMatch> matches = HaplotyperTest.testCallHaplotype(m_tsvFile, vcfFile, true, false, true);

    List<String> expectedMatches = Lists.newArrayList("*1A/*15","*1B/*5");
    TestUtil.assertDiplotypePairs(expectedMatches, matches);
  }

  @Test
  public void slco1b1s1as15s1bs5TwoMissing() throws Exception {
    /* Test two positions missing:
      chr12	21239145	rs200995543	C	T	.	PASS	only--star-3-4	GT	0/0
      chr12	21239158	rs140790673	C	T	.	PASS	second-star-29	GT	0/0


    *1a/*15. Except we can't distinguish *1B/*5. *5/*29 Matches as well,
     because the first position matches.  However this shows difference from
     single missing position. In the output html we see:

       There were 2 missing positions from the VCF file:

        21239158 (g.21239158C>T)
        21239145 (g.21239145C>T)
        The following haplotype(s) were eliminated from consideration:

        *34

    This shows difference in behaviour  between *34 (a single snp) and  *29 (two snps). *29 is not excluded. Again need to
    discus with others to confirm this logic.

     */

    Path vcfFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/SLCO1B1/s1as15s1bs5twomissing.vcf");
    List<DiplotypeMatch> matches = HaplotyperTest.testCallHaplotype(m_tsvFile, vcfFile, true, false, true);

    List<String> expectedMatches = Lists.newArrayList("*1A/*15","*1B/*5");
    TestUtil.assertDiplotypePairs(expectedMatches, matches);
  }


  @Test
  public void s1as15s1bs5NotCalled() throws Exception {
    // Test *1a/*15. Except we can't distinguish *1B/*5 or *5/*29, because final position is ./.

    Path vcfFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/SLCO1B1/s1as15s1bs5notcalled.vcf");
    List<DiplotypeMatch> matches = HaplotyperTest.testCallHaplotype(m_tsvFile, vcfFile, true, false, true);

    // same as above, but using ./. to signify missing snp.
    List<String> expectedMatches = Lists.newArrayList("*1A/*15","*1B/*5");
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

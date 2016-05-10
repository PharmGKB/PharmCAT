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
public class HaplotyperDpydTest {
  private Path m_jsonFile;

  @Before
  public void before() throws Exception {
    m_jsonFile =  TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/DPYD_translation.json");
  }


  @Test
  public void dpyds1s1() throws Exception {
    // Test *1/*1 - contains a del
    /* TODO: Lester/Mark - Fails. In the output vcf for *1*1 there is no call for rs72549303 position ("|")
    This issue is also seen in cftr F508d/F508d example.  Here there is two  dels

    Also, per my understanding, we are not looking for ATGA in *1 for rs72549309, but just A

    */

    Path vcfFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/DPYD/s1s1.vcf");
    List<DiplotypeMatch> matches = HaplotyperTest.testCallHaplotype(m_jsonFile, vcfFile);

    List<String> expectedMatches = Lists.newArrayList("*1/*1");
    TestUtil.assertDiplotypePairs(expectedMatches, matches);
  }


  @Test
  public void dpyds2aRs67376798A() throws Exception {
    // Test s2aRs67376798A
    /* TODO: Lester/Mark - Fails.

    There is a '/' in the star name. Could cause issues?
    Expected: [*2A/rs67376798T/A]
    Got:      []

    Again no vcf calls in rs72549303 or 	rs72549309 so probably failing because of that

    */

    Path vcfFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/DPYD/s2aRs67376798A.vcf");
    List<DiplotypeMatch> matches = HaplotyperTest.testCallHaplotype(m_jsonFile, vcfFile);

    List<String> expectedMatches = Lists.newArrayList("*2A/rs67376798T/A");
    TestUtil.assertDiplotypePairs(expectedMatches, matches);
  }

  @Test
  public void dpyds1s7() throws Exception {
    // Test *1/*7
    /* TODO: Lester/Mark - Fails. Combination of above issue (blanks in vcf match) and below issue:

    IMPORTANT! - this is how I think the vcf should be represented:
    chr1	97740414	rs72549308	GATGA	G	.	PASS	assume-default	GT	0/1
    chr1	97740415	rs72549309	A	.	.	PASS	assume-default	GT	0/0

    Notice the added 97740414 position for the deletion, which is not in the tsv. So -1 rule
    for deletions, with the -1 nucleotide anchoring the deletion. This is based on the 4.2 spec "For
    simple insertions and deletions in which either the REF or one of the ALT alleles would otherwise
    be null/empty, the REF and ALT Strings must include the base before the event (which must be
    reflected in the POS field)".

    My reading - for deletions we need to look at position -1. For insertions the -1 chr to where the insertions start
    is the equivalent to where the rsid is already, so don't need to go back a position.

    */

    Path vcfFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/DPYD/s1s7.vcf");
    List<DiplotypeMatch> matches = HaplotyperTest.testCallHaplotype(m_jsonFile, vcfFile);

    List<String> expectedMatches = Lists.newArrayList("*1/*7");
    TestUtil.assertDiplotypePairs(expectedMatches, matches);
  }

}

package org.pharmgkb.pharmcat.haplotype;

import java.nio.file.Path;
import java.util.List;
import com.google.common.collect.Lists;
import org.junit.Before;
import org.junit.Test;
import org.pharmgkb.common.util.PathUtils;
import org.pharmgkb.pharmcat.haplotype.model.Result;

import static org.pharmgkb.pharmcat.haplotype.NamedAlleleMatcherTest.assertDiplotypePairs;
import static org.pharmgkb.pharmcat.haplotype.NamedAlleleMatcherTest.testMatchNamedAlleles;


/**
 * JUnit test for {@link NamedAlleleMatcher#callDiplotypes(MatchData)}.
 *
 * @author Lester Carter
 */
public class NamedAlleleMatcherDpydTest {
  private Path m_definitionFile;

  @Before
  public void before() throws Exception {
    m_definitionFile =  PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/DPYD_translation.json");
  }


  @Test
  public void dpyds1s1() throws Exception {
    // Test *1/*1 - contains a del

    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/DPYD/s1s1.vcf");
    List<String> expectedMatches = Lists.newArrayList("*1/*1");

    Result result = testMatchNamedAlleles(m_definitionFile, vcfFile);
    assertDiplotypePairs(expectedMatches, result);
  }


  @Test
  public void dpyds2aRs67376798A() throws Exception {
    // Test *2a/Rs67376798A

    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/DPYD/s2aRs67376798A.vcf");
    List<String> expectedMatches = Lists.newArrayList("*2A/rs67376798T/A");

    Result result = testMatchNamedAlleles(m_definitionFile, vcfFile);
    assertDiplotypePairs(expectedMatches, result);
  }

  @Test
  public void dpyds1s2b() throws Exception {
    // Test *1/*2b - however can't be distinguished from *2A/*5

    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/DPYD/s1s2b.vcf");
    List<String> expectedMatches = Lists.newArrayList("*1/*2B", "*2A/*5");

    Result result = testMatchNamedAlleles(m_definitionFile, vcfFile);
    assertDiplotypePairs(expectedMatches, result);
  }

  @Test
  public void dpyds1s7() throws Exception {
    // Test *1/*7

    /*
    IMPORTANT TEST CASE! - this is how I think the vcf should be represented:
    chr1	97740414	rs72549308	GATGA	G	.	PASS	assume-default	GT	0/1
    chr1	97740415	rs72549309	A	.	.	PASS	assume-default	GT	0/0

    Notice the added 97740414 position for the deletion, which is not in the tsv. So -1 rule
    for deletions, with the -1 nucleotide anchoring the deletion. This is based on the 4.2 spec "For
    simple insertions and deletions in which either the REF or one of the ALT alleles would otherwise
    be null/empty, the REF and ALT Strings must include the base before the event (which must be
    reflected in the POS field)".

    My reading - for deletions we need to look at position -1. For insertions the -1 position to where the insertions start
    is the equivalent to where the rsid is already, so don't need to go back a position.

    */

    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/DPYD/s1s7.vcf");
    List<String> expectedMatches = Lists.newArrayList("*1/*7");

    Result result = testMatchNamedAlleles(m_definitionFile, vcfFile);
    assertDiplotypePairs(expectedMatches, result);
  }
}

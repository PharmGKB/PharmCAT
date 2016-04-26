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
 * @author Mark Woon
 */
public class HaplotyperCyp2c19Test {
  private Path m_tsvFile;

  @Before
  public void before() throws Exception {
    m_tsvFile =  TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/CYP2C19.tsv");
  }


  private List<DiplotypeMatch> testCallHaplotype(Path vcfFile) throws Exception {

    DefinitionReader definitionReader = new DefinitionReader();
    definitionReader.read(m_tsvFile);
    String gene = definitionReader.getHaplotypes().keySet().iterator().next();

    VcfReader vcfReader = new VcfReader(Haplotyper.calculateLocationsOfInterest(definitionReader));
    SortedMap<String, SampleAllele> alleleMap = vcfReader.read(vcfFile);

    Haplotyper haplotyper = new Haplotyper(definitionReader);
    List<DiplotypeMatch> matches = haplotyper.callDiplotypes(alleleMap, gene);
    StringBuilder rezBuilder = new StringBuilder();
    for (DiplotypeMatch dm : matches) {
      if (rezBuilder.length() > 0) {
        rezBuilder.append(", ");
      }
      rezBuilder.append(dm.getName())
          .append(" (")
          .append(dm.getScore())
          .append(")");
    }
    System.out.println(rezBuilder);

    // print
    new Report(definitionReader)
        .forFile(vcfFile)
        .gene(gene, matches, alleleMap.values())
        .printHtml();

    return matches;
  }


  @Test
  public void cyp2c19s1s1() throws Exception {
    // Simple *1/*1 reference test
    Path vcfFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/cyp2c19/s1s1.vcf");
    List<DiplotypeMatch> matches = testCallHaplotype(vcfFile);

    List<String> expectedMatches = Lists.newArrayList("*1/*1");
    TestUtil.assertDiplotypePairs(expectedMatches, matches);
  }


  @Test
  public void cyp2c19s1s2() throws Exception {
    // Test simple case of one heterozygous snp
    Path vcfFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/cyp2c19/s1s2.vcf");
    List<DiplotypeMatch> matches = testCallHaplotype(vcfFile);

    List<String> expectedMatches = Lists.newArrayList("*1/*2");
    TestUtil.assertDiplotypePairs(expectedMatches, matches);
  }


  @Test
  public void cyp2c19s1s4b() throws Exception {
    /* Test *1 *4b. s1s4b-longest wins(not s4as17)
    // TODO: still in progress
    Expected: [*1/*4b]
    Got:      [*1/*17, *1/*4A, *1/*4B, *4A/*17]

    *1/*17 is impossible - if *17 then rs12248560 and rs28399504 must be in trans, so *1 is excluded, and it must
    be 4a/*17.

    If rs12248560 and rs28399504 are in cis then *4B and other copy is wild type, hence *1/*4b

    The options are [*1/*4B, *4A/*17], but *1/*4B is the longest haplotype match (2 matches)

    Michelle would expect *1/4B

    */
    Path vcfFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/cyp2c19/s1s4b.vcf");
    List<DiplotypeMatch> matches = testCallHaplotype(vcfFile);

    List<String> expectedMatches = Lists.newArrayList("*1/*4B");
    TestUtil.assertDiplotypePairs(expectedMatches, matches);
  }


  @Test
  public void cyp2c19s1s28() throws Exception {
    // Test *1 *28. The longest possible match should win.  *15 is possible but should be excluded by longest match
    // *1/*15 would be excluded because if *15 the other snps would preclude *1.
    // TODO: still in progress
    Path vcfFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/cyp2c19/s1s28.vcf");
    List<DiplotypeMatch> matches = testCallHaplotype(vcfFile);

    List<String> expectedMatches = Lists.newArrayList("*1/*28");
    TestUtil.assertDiplotypePairs(expectedMatches, matches);
  }


  @Test
  public void cyp2c19s2s2() throws Exception {
    // Test simple case of one homozygous snp
    Path vcfFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/cyp2c19/s2s2.vcf");
    List<DiplotypeMatch> matches = testCallHaplotype(vcfFile);

    List<String> expectedMatches = Lists.newArrayList("*2/*2");
    TestUtil.assertDiplotypePairs(expectedMatches, matches);
  }


  @Test
  public void cyp2c19s2s3() throws Exception {
    /* Test case of two snps - *2 and *3
    TODO: still in progress
    Initial *2/*2 bug fixed.  However:

    Expected: [*2/*3]
    Got:      [*1/*2, *1/*3, *2/*3]

    If *2 is on one copy, then the *3 must be on the other, excluding *1.
     */
    Path vcfFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/cyp2c19/s2s3.vcf");
    List<DiplotypeMatch> matches = testCallHaplotype(vcfFile);

    List<String> expectedMatches = Lists.newArrayList("*2/*3");
    TestUtil.assertDiplotypePairs(expectedMatches, matches);
  }


  @Test
  public void cyp2c19s4as4b() throws Exception {
    /* Test *4a *4b. het and homo rsids
    // TODO: still in progress
    Expected: [*4a/*4b]
    Got:      [*4A/*17, *4A/*4A, *4A/*4B]

    rs12248560 (*4b and *17) is het
    rs28399504 (4b and 4a) is homo

    Possible options;
          rs12248560 (C/T)  rs28399504 (G/G)
               1 2              1 2
               C T              G G     *4a *4b
               T C              G G     *4b *4a (same as above)
    Can't be any of the following:
     if *4A*17 C  T              G G     *4a - G in first position means T must be second, *17 matches, but *4b is longest
     IF *4a*4a C !T              G G     *4a - G in first position means T must be second, second *2
    */
    Path vcfFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/cyp2c19/s4as4b.vcf");
    List<DiplotypeMatch> matches = testCallHaplotype(vcfFile);

    List<String> expectedMatches = Lists.newArrayList("*4A/*4B");
    TestUtil.assertDiplotypePairs(expectedMatches, matches);
  }


  @Test
  public void cyp2c19s4bs17() throws Exception {
    /* Test *4b/*17
    TODO: still in progress
    Expected: [*4b/*17]
    Got:      [*17/*17, *4A/*17, *4B/*17]

    *17/*17 diplotype not possible due to rs28399504 het call
    * 4A/*17 not possible because rs12248560 is homo, so 4a will have to become the 4b call

    */
    Path vcfFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/cyp2c19/s4bs17.vcf");
    List<DiplotypeMatch> matches = testCallHaplotype(vcfFile);

    List<String> expectedMatches = Lists.newArrayList("*4B/*17");
    TestUtil.assertDiplotypePairs(expectedMatches, matches);
  }


  @Test
  public void cyp2c19s15s28() throws Exception {
    // Test *15 *28. The shared position is homo
    /* TODO: still in progress

    Expected: [*15/*28]
    Got:      [*15/*15, *15/*28]

    *15/*15 is wrong here - if *15 then all the *28 positions would mean it can't be another *15.

     */

    Path vcfFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/cyp2c19/s15s28.vcf");
    List<DiplotypeMatch> matches = testCallHaplotype(vcfFile);

    List<String> expectedMatches = Lists.newArrayList("*15/*28");
    TestUtil.assertDiplotypePairs(expectedMatches, matches);
  }


  //@Test
  public void cyp2c19sUnks17() throws Exception {
    /* Test *Unk/*17 - This is the case when one haplotype matches, but the second would be undefined
    TODO: still in progress
    Expected: [*Unk/*17]
    Got:      [*17/*17]

    Can't be *17/*17
          rs12248560 (T/T)  rs28399504 (A/C)
            1   2               1   2
            T   T               A   C   TA matches *17, however TC matches nothing in table
            T   T               C   A   As above
    */

    Path vcfFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/cyp2c19/sUnks17.vcf");
    List<DiplotypeMatch> matches = testCallHaplotype(vcfFile);

    List<String> expectedMatches = Lists.newArrayList("*Unk/*1","*Unk/*17");
    TestUtil.assertDiplotypePairs(expectedMatches, matches);
  }
}

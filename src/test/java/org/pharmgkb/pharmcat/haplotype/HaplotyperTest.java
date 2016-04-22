package org.pharmgkb.pharmcat.haplotype;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Collection;
import java.util.List;
import java.util.Set;
import java.util.SortedMap;
import com.google.common.collect.Sets;
import org.junit.Test;
import org.pharmgkb.pharmcat.TestUtil;
import org.pharmgkb.pharmcat.haplotype.model.DiplotypeMatch;


/**
 * JUnit test for {@link Haplotyper}.
 *
 * @author Mark Woon
 */
public class HaplotyperTest {

  @Test
  public void testCyp2c19() throws Exception {
    // initial test - file reading etc
    Path tsvFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/CYP2C19.tsv");
    Path vcfFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/NA12878.2c19_filtered.vcf");

    DefinitionReader definitionReader = new DefinitionReader();
    definitionReader.read(tsvFile);

    Haplotyper haplotyper = new Haplotyper(definitionReader);
    Path jsonFile = Files.createTempFile("haplotyper", ".json");
    Report report = haplotyper.call(vcfFile);
  }


  @Test
  public void cyp2c19s1s1() throws Exception {
    // Simple *1/*1 reference test
    Path tsvFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/CYP2C19.tsv");
    Path vcfFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/HaplotyperTest/NA12878.cyp2c19.s1s1.vcf");
    String gene = "CYP2C19";

    DefinitionReader definitionReader = new DefinitionReader();
    definitionReader.read(tsvFile);

    VcfReader vcfReader = new VcfReader(Haplotyper.calculateLocationsOfInterest(definitionReader));

    Haplotyper haplotyper = new Haplotyper(definitionReader);
    SortedMap<String, SampleAllele> alleles = vcfReader.read(vcfFile);

    List<DiplotypeMatch> matches = haplotyper.callHaplotype(alleles, gene);
    Set<String> expectedMatches = Sets.newHashSet("*1/*1");
    TestUtil.assertDiplotypePairs(expectedMatches, matches);
    System.out.println(matches);
    ComparisonUtil.printMatchPairs(matches);
  }



  @Test
  public void cyp2c19s1s2() throws Exception {
    // Test simple case of one heterozygous snp
    Path tsvFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/CYP2C19.tsv");
    Path vcfFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/HaplotyperTest/HaplotyperTest.cyp2c19s1s2.vcf");
    String gene = "CYP2C19";

    DefinitionReader definitionReader = new DefinitionReader();
    definitionReader.read(tsvFile);

    VcfReader vcfReader = new VcfReader(Haplotyper.calculateLocationsOfInterest(definitionReader));

    Haplotyper haplotyper = new Haplotyper(definitionReader);
    SortedMap<String, SampleAllele> alleles = vcfReader.read(vcfFile);

    List<DiplotypeMatch> matches = haplotyper.callHaplotype(alleles, gene);
    Set<String> expectedMatches = Sets.newHashSet("*1/*2");
    TestUtil.assertDiplotypePairs(expectedMatches, matches);
    System.out.println(matches);
  }

  @Test
  public void cyp2c19s2s2() throws Exception {
    // Test simple case of one homozygous snp
    Path tsvFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/CYP2C19.tsv");
    Path vcfFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/HaplotyperTest/HaplotyperTest.cyp2c19s2s2.vcf");
    String gene = "CYP2C19";

    DefinitionReader definitionReader = new DefinitionReader();
    definitionReader.read(tsvFile);

    VcfReader vcfReader = new VcfReader(Haplotyper.calculateLocationsOfInterest(definitionReader));

    Haplotyper haplotyper = new Haplotyper(definitionReader);
    SortedMap<String, SampleAllele> alleles = vcfReader.read(vcfFile);

    List<DiplotypeMatch> matches = haplotyper.callHaplotype(alleles, gene);
    Set<String> expectedMatches = Sets.newHashSet("*2/*2");
    TestUtil.assertDiplotypePairs(expectedMatches, matches);
    System.out.println(matches);
    ComparisonUtil.printMatchPairs(matches);
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
    Path tsvFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/CYP2C19.tsv");
    Path vcfFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/HaplotyperTest/HaplotyperTest.cyp2c19s2s3.vcf");
    String gene = "CYP2C19";

    DefinitionReader definitionReader = new DefinitionReader();
    definitionReader.read(tsvFile);

    VcfReader vcfReader = new VcfReader(Haplotyper.calculateLocationsOfInterest(definitionReader));

    Haplotyper haplotyper = new Haplotyper(definitionReader);
    SortedMap<String, SampleAllele> alleleMap = vcfReader.read(vcfFile);

    List<DiplotypeMatch> matches = haplotyper.callHaplotype(alleleMap, gene);
    printReport(definitionReader, vcfFile, gene, matches, alleleMap.values());

    Set<String> expectedMatches = Sets.newHashSet("*2/*3");
    TestUtil.assertDiplotypePairs(expectedMatches, matches);
    System.out.println(matches);
  }


  @Test
  public void cyp2c19s1s28() throws Exception {
    // Test *1 *28. The longest possible match should win.  *15 is possible but should be excluded by longest match
    // *1/*15 would be excluded because if *15 the other snps would preclude *1.
    // TODO: still in progress
    Path tsvFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/CYP2C19.tsv");
    Path vcfFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/HaplotyperTest/HaplotyperTest.cyp2c19s1s28.vcf");
    String gene = "CYP2C19";

    DefinitionReader definitionReader = new DefinitionReader();
    definitionReader.read(tsvFile);

    VcfReader vcfReader = new VcfReader(Haplotyper.calculateLocationsOfInterest(definitionReader));

    Haplotyper haplotyper = new Haplotyper(definitionReader);
    SortedMap<String, SampleAllele> alleles = vcfReader.read(vcfFile);

    List<DiplotypeMatch> matches = haplotyper.callHaplotype(alleles, gene);
    Set<String> expectedMatches = Sets.newHashSet("*1/*28");
    TestUtil.assertDiplotypePairs(expectedMatches, matches);
    System.out.println(matches);
  }


  @Test
  public void cyp2c19s15s28() throws Exception {
    // Test *15 *28. The shared position is homo
    /* TODO: still in progress

Expected: [*15/*28]
Got:      [*15/*15, *15/*28]

    *15/*15 is wrong here - if *15 then all the *28 positions would mean it can't be another *15.

     */

    Path tsvFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/CYP2C19.tsv");
    Path vcfFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/HaplotyperTest/HaplotyperTest.cyp2c19s15s28.vcf");
    String gene = "CYP2C19";

    DefinitionReader definitionReader = new DefinitionReader();
    definitionReader.read(tsvFile);

    VcfReader vcfReader = new VcfReader(Haplotyper.calculateLocationsOfInterest(definitionReader));

    Haplotyper haplotyper = new Haplotyper(definitionReader);
    SortedMap<String, SampleAllele> alleles = vcfReader.read(vcfFile);

    List<DiplotypeMatch> matches = haplotyper.callHaplotype(alleles, gene);
    Set<String> expectedMatches = Sets.newHashSet("*15/*28");
    TestUtil.assertDiplotypePairs(expectedMatches, matches);
    System.out.println(matches);
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

    Path tsvFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/CYP2C19.tsv");
    Path vcfFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/HaplotyperTest/HaplotyperTest.cyp2c19s1s4b.vcf");
    String gene = "CYP2C19";

    DefinitionReader definitionReader = new DefinitionReader();
    definitionReader.read(tsvFile);

    VcfReader vcfReader = new VcfReader(Haplotyper.calculateLocationsOfInterest(definitionReader));

    Haplotyper haplotyper = new Haplotyper(definitionReader);
    SortedMap<String, SampleAllele> alleles = vcfReader.read(vcfFile);

    List<DiplotypeMatch> matches = haplotyper.callHaplotype(alleles, gene);
    Set<String> expectedMatches = Sets.newHashSet("*1/*4b");
    TestUtil.assertDiplotypePairs(expectedMatches, matches);
    System.out.println(matches);
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

    Path tsvFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/CYP2C19.tsv");
    Path vcfFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/HaplotyperTest/HaplotyperTest.cyp2c19s4as4b.vcf");
    String gene = "CYP2C19";

    DefinitionReader definitionReader = new DefinitionReader();
    definitionReader.read(tsvFile);

    VcfReader vcfReader = new VcfReader(Haplotyper.calculateLocationsOfInterest(definitionReader));

    Haplotyper haplotyper = new Haplotyper(definitionReader);
    SortedMap<String, SampleAllele> alleles = vcfReader.read(vcfFile);

    List<DiplotypeMatch> matches = haplotyper.callHaplotype(alleles, gene);
    Set<String> expectedMatches = Sets.newHashSet("*4a/*4b");
    TestUtil.assertDiplotypePairs(expectedMatches, matches);
    System.out.println(matches);
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

    Path tsvFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/CYP2C19.tsv");
    Path vcfFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/HaplotyperTest/HaplotyperTest.cyp2c19s4bs17.vcf");
    String gene = "CYP2C19";

    DefinitionReader definitionReader = new DefinitionReader();
    definitionReader.read(tsvFile);

    VcfReader vcfReader = new VcfReader(Haplotyper.calculateLocationsOfInterest(definitionReader));

    Haplotyper haplotyper = new Haplotyper(definitionReader);
    SortedMap<String, SampleAllele> alleles = vcfReader.read(vcfFile);

    List<DiplotypeMatch> matches = haplotyper.callHaplotype(alleles, gene);
    Set<String> expectedMatches = Sets.newHashSet("*4b/*17");
    TestUtil.assertDiplotypePairs(expectedMatches, matches);
    System.out.println(matches);
  }


  @Test
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

    Path tsvFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/CYP2C19.tsv");
    Path vcfFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/HaplotyperTest/HaplotyperTest.cyp2c19sUnks17.vcf");
    String gene = "CYP2C19";

    DefinitionReader definitionReader = new DefinitionReader();
    definitionReader.read(tsvFile);

    VcfReader vcfReader = new VcfReader(Haplotyper.calculateLocationsOfInterest(definitionReader));

    Haplotyper haplotyper = new Haplotyper(definitionReader);
    SortedMap<String, SampleAllele> alleles = vcfReader.read(vcfFile);

    List<DiplotypeMatch> matches = haplotyper.callHaplotype(alleles, gene);
    Set<String> expectedMatches = Sets.newHashSet("*Unk/*1","*Unk/*17");
    TestUtil.assertDiplotypePairs(expectedMatches, matches);
    System.out.println(matches);
  }



  private void printReport(DefinitionReader definitionReader, Path vcfFile, String gene,
      List<DiplotypeMatch> matches, Collection<SampleAllele> alleles) throws IOException {

    new Report(definitionReader)
        .forFile(vcfFile)
        .gene(gene, matches, alleles)
        .printHtml();
  }

}

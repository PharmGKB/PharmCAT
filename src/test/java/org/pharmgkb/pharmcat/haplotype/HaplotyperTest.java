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
    haplotyper.call(vcfFile, jsonFile);
  }


  @Test
  public void cyp2c19s1s1() throws Exception {
    // Simple *1/*1 reference test
    Path tsvFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/CYP2C19.tsv");
    Path vcfFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/NA12878.cyp2c19.s1s1.vcf");
    String gene = "CYP2C19";

    DefinitionReader definitionReader = new DefinitionReader();
    definitionReader.read(tsvFile);

    VcfReader vcfReader = new VcfReader(Haplotyper.calculateLocationsOfInterest(definitionReader));

    Haplotyper haplotyper = new Haplotyper(definitionReader);
    SortedMap<String, SampleAllele> alleles = vcfReader.read(vcfFile);

    List<List<HaplotypeMatch>> matches = haplotyper.callHaplotype(alleles, gene);
    Set<String> expectedMatches = Sets.newHashSet("*1/*1");
    TestUtil.assertDiplotypePairs(expectedMatches, matches);
    System.out.println(matches);
    ComparisonUtil.printMatchPairs(matches);
  }



  @Test
  public void cyp2c19s1s2() throws Exception {
    // Test simple case of one heterozygous snp
    Path tsvFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/CYP2C19.tsv");
    Path vcfFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/HaplotyperTest.cyp2c19s1s2.vcf");
    String gene = "CYP2C19";

    DefinitionReader definitionReader = new DefinitionReader();
    definitionReader.read(tsvFile);

    VcfReader vcfReader = new VcfReader(Haplotyper.calculateLocationsOfInterest(definitionReader));

    Haplotyper haplotyper = new Haplotyper(definitionReader);
    SortedMap<String, SampleAllele> alleles = vcfReader.read(vcfFile);

    List<List<HaplotypeMatch>> matches = haplotyper.callHaplotype(alleles, gene);
    Set<String> expectedMatches = Sets.newHashSet("*1/*2");
    TestUtil.assertDiplotypePairs(expectedMatches, matches);
    System.out.println(matches);
    ComparisonUtil.printMatchPairs(matches);
  }

  @Test
  public void cyp2c19s1s2Homo() throws Exception {
    // Test simple case of one homozygous snp
    Path tsvFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/CYP2C19.tsv");
    Path vcfFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/HaplotyperTest.cyp2c19s1s2homo.vcf");
    String gene = "CYP2C19";

    DefinitionReader definitionReader = new DefinitionReader();
    definitionReader.read(tsvFile);

    VcfReader vcfReader = new VcfReader(Haplotyper.calculateLocationsOfInterest(definitionReader));

    Haplotyper haplotyper = new Haplotyper(definitionReader);
    SortedMap<String, SampleAllele> alleles = vcfReader.read(vcfFile);

    List<List<HaplotypeMatch>> matches = haplotyper.callHaplotype(alleles, gene);
    Set<String> expectedMatches = Sets.newHashSet("*2/*2");
    TestUtil.assertDiplotypePairs(expectedMatches, matches);
    System.out.println(matches);
    ComparisonUtil.printMatchPairs(matches);
  }


  @Test
  public void cyp2c19s2s3() throws Exception {
    /* Test case of two snps - *2 and *3
    TODO:
    Below makes no sense - can't be *2/*2 etc when het

    Expected: [*2/*3]
    Got:      [*1/*2, *1/*3, *2/*2, *2/*3, *3/*3]

    If *3 is one one allele the *2 has to be on the other. If there were both on the same then it would be another allele?  Presuming that *1/*2+*3 is not valid
    It's possible that make is using *1 in place of *unk?

     */

    Path tsvFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/CYP2C19.tsv");
    Path vcfFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/HaplotyperTest.cyp2c19s2s3.vcf");
    String gene = "CYP2C19";

    DefinitionReader definitionReader = new DefinitionReader();
    definitionReader.read(tsvFile);

    VcfReader vcfReader = new VcfReader(Haplotyper.calculateLocationsOfInterest(definitionReader));

    Haplotyper haplotyper = new Haplotyper(definitionReader);
    SortedMap<String, SampleAllele> alleles = vcfReader.read(vcfFile);

    List<List<HaplotypeMatch>> matches = haplotyper.callHaplotype(alleles, gene);
    printReport(definitionReader, vcfFile, gene, matches, alleles.values());
    Set<String> expectedMatches = Sets.newHashSet("*2/*3");
    TestUtil.assertDiplotypePairs(expectedMatches, matches);
    System.out.println(matches);
    ComparisonUtil.printMatchPairs(matches);
  }

  @Test
  public void cyp2c19s15s28() throws Exception {
    // Test *15 *28. The shared position is het so it's either *1/*28(?) or ["*1/*28","*unk, *15"] - not MW's preference
    // TODO: // still work in progress
    Path tsvFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/CYP2C19.tsv");
    Path vcfFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/HaplotyperTest.cyp2c19s15s28.vcf");
    String gene = "CYP2C19";

    DefinitionReader definitionReader = new DefinitionReader();
    definitionReader.read(tsvFile);

    VcfReader vcfReader = new VcfReader(Haplotyper.calculateLocationsOfInterest(definitionReader));

    Haplotyper haplotyper = new Haplotyper(definitionReader);
    SortedMap<String, SampleAllele> alleles = vcfReader.read(vcfFile);

    List<List<HaplotypeMatch>> matches = haplotyper.callHaplotype(alleles, gene);
    Set<String> expectedMatches = Sets.newHashSet("*4b/*17");
    //TestUtil.assertDiplotypePairs(expectedMatches, matches);
    System.out.println(matches);
    ComparisonUtil.printMatchPairs(matches);
  }


  @Test
  public void cyp2c19s4bs17() throws Exception {
    // Test differentiation between *4b/*17
    // TODO: still working in progress
    Path tsvFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/CYP2C19.tsv");
    Path vcfFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/HaplotyperTest.cyp2c19s4bs17.vcf");
    String gene = "CYP2C19";

    DefinitionReader definitionReader = new DefinitionReader();
    definitionReader.read(tsvFile);

    VcfReader vcfReader = new VcfReader(Haplotyper.calculateLocationsOfInterest(definitionReader));

    Haplotyper haplotyper = new Haplotyper(definitionReader);
    SortedMap<String, SampleAllele> alleles = vcfReader.read(vcfFile);

    List<List<HaplotypeMatch>> matches = haplotyper.callHaplotype(alleles, gene);
    Set<String> expectedMatches = Sets.newHashSet("*4b/*17");
    // TestUtil.assertDiplotypePairs(expectedMatches, matches);
    // System.out.println(matches);
    // ComparisonUtil.printMatchPairs(matches);
  }


  private void printReport(DefinitionReader definitionReader, Path vcfFile, String gene,
      List<List<HaplotypeMatch>> matches, Collection<SampleAllele> alleles) throws IOException {

    new JsonReport(definitionReader)
        .forFile(vcfFile)
        .haplotype(gene, matches, alleles)
        .printHtml();
  }

}

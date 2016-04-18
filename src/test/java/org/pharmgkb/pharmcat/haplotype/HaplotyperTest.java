package org.pharmgkb.pharmcat.haplotype;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Map;
import java.util.Set;
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

    Path tsvFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/CYP2C19.tsv");
    Path vcfFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/NA12878.2c19_filtered.vcf");

    DefinitionReader definitionReader = new DefinitionReader();
    definitionReader.read(tsvFile);

    Haplotyper haplotyper = new Haplotyper(definitionReader);
    Path jsonFile = Files.createTempFile("haplotyper", ".json");
    haplotyper.call(vcfFile, jsonFile);
  }


  @Test
  public void test1() throws Exception {

    Path tsvFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/CYP2C19.tsv");
    Path vcfFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/NA12878.2c19_filtered.vcf");
    String gene = "CYP2C19";

    DefinitionReader definitionReader = new DefinitionReader();
    definitionReader.read(tsvFile);

    VcfReader vcfReader = new VcfReader(Haplotyper.calculateLocationsOfInterest(definitionReader));

    Haplotyper haplotyper = new Haplotyper(definitionReader);
    Map<String, SampleAllele> alleles = vcfReader.read(vcfFile);

    List<List<HaplotypeMatch>> matches = haplotyper.callHaplotype(alleles, gene);
    Set<String> expectedMatches = Sets.newHashSet("*1/*3", "*1/*27", "*3/*3", "*3/*27", "*27/*27");
    TestUtil.assertDiplotypePairs(expectedMatches, matches);
    System.out.println(matches);
    ComparisonUtil.printMatchPairs(matches);
  }


  @Test
  public void cyp2c19s1s1() throws Exception {
    // Fails no matches. WARN VcfReader - Unexpected state: no alt alleles for chr10:94760676 etc
    Path tsvFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/CYP2C19.tsv");
    Path vcfFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/NA12878.cyp2c19.s1s1.vcf");
    String gene = "CYP2C19";

    DefinitionReader definitionReader = new DefinitionReader();
    definitionReader.read(tsvFile);

    VcfReader vcfReader = new VcfReader(Haplotyper.calculateLocationsOfInterest(definitionReader));

    Haplotyper haplotyper = new Haplotyper(definitionReader);
    Map<String, SampleAllele> alleles = vcfReader.read(vcfFile);

    List<List<HaplotypeMatch>> matches = haplotyper.callHaplotype(alleles, gene);
    Set<String> expectedMatches = Sets.newHashSet("*1/*1");
    TestUtil.assertDiplotypePairs(expectedMatches, matches);
    System.out.println(matches);
    ComparisonUtil.printMatchPairs(matches);
  }

  @Test
  public void cyp2c19s1s1Alt() throws Exception {
    // Fails with no matches
    Path tsvFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/CYP2C19.tsv");
    Path vcfFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/NA12878.cyp2c19.s1s1.alt.vcf");
    String gene = "CYP2C19";

    DefinitionReader definitionReader = new DefinitionReader();
    definitionReader.read(tsvFile);

    VcfReader vcfReader = new VcfReader(Haplotyper.calculateLocationsOfInterest(definitionReader));

    Haplotyper haplotyper = new Haplotyper(definitionReader);
    Map<String, SampleAllele> alleles = vcfReader.read(vcfFile);

    List<List<HaplotypeMatch>> matches = haplotyper.callHaplotype(alleles, gene);
    Set<String> expectedMatches = Sets.newHashSet("*1/*1");
    TestUtil.assertDiplotypePairs(expectedMatches, matches);
    System.out.println(matches);
    ComparisonUtil.printMatchPairs(matches);
  }

}

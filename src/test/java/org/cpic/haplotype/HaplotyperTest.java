package org.cpic.haplotype;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Map;
import org.cpic.TestUtil;
import org.junit.Test;


/**
 * JUnit test for {@link Haplotyper}.
 *
 * @author Mark Woon
 */
public class HaplotyperTest {

  @Test
  public void testCyp2c19() throws Exception {

    Path tsvFile = TestUtil.getFile("org/cpic/haplotype/CYP2C19.tsv");
    Path vcfFile = TestUtil.getFile("org/cpic/haplotype/NA12878.2c19_filtered.vcf");

    DefinitionReader definitionReader = new DefinitionReader();
    definitionReader.read(tsvFile);

    Haplotyper haplotyper = new Haplotyper(definitionReader);
    Path jsonFile = Files.createTempFile("haplotyper", ".json");
    haplotyper.call(vcfFile, jsonFile);
  }


  @Test
  public void test1() throws Exception {

    Path tsvFile = TestUtil.getFile("org/cpic/haplotype/CYP2C19.tsv");
    Path vcfFile = TestUtil.getFile("org/cpic/haplotype/NA12878.2c19_filtered.vcf");
    String gene = "CYP2C19";

    DefinitionReader definitionReader = new DefinitionReader();
    definitionReader.read(tsvFile);

    VcfReader vcfReader = new VcfReader(Haplotyper.calculateLocationsOfInterest(definitionReader));

    Haplotyper haplotyper = new Haplotyper(definitionReader);
    for (String key : definitionReader.getHaplotypes().keySet()) {
      for (Haplotype hap : definitionReader.getHaplotypes().get(key)) {
        System.out.println(hap.getPermutations().pattern());
      }
    }
    Map<String, SampleAllele> alleles = vcfReader.read(vcfFile);

    List<List<HaplotypeMatch>> matches = haplotyper.callHaplotype(alleles, gene);
    System.out.println(matches);
  }
}

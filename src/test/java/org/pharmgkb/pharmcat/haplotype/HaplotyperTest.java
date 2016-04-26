package org.pharmgkb.pharmcat.haplotype;

import java.nio.file.Path;
import org.junit.Test;
import org.pharmgkb.pharmcat.TestUtil;


/**
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
    Report report = haplotyper.call(vcfFile);
  }
}

package org.pharmgkb.pharmcat.haplotype;

import java.nio.file.Path;
import org.junit.jupiter.api.Test;
import org.pharmgkb.pharmcat.TestVcfBuilder;
import org.pharmgkb.pharmcat.util.DataManager;

import static org.pharmgkb.pharmcat.haplotype.NamedAlleleMatcherTest.assertDiplotypePairs;
import static org.pharmgkb.pharmcat.haplotype.NamedAlleleMatcherTest.testMatchNamedAlleles;


/**
 * Test calling for CYP4F2.
 *
 * @author Mark Woon
 */
class NamedAlleleMatcherCyp4f2Test {
  private final Path sf_definitionFile = DataManager.DEFAULT_DEFINITION_DIR.resolve("CYP4F2_translation.json");


  @Test
  void s2_s3() throws Exception {
    Path vcfFile = new TestVcfBuilder("*2/*3")
        .reference("CYP4F2")
        .saveFile()
        .variation("CYP4F2", "rs3093105", "C", "A")
        .variation("CYP4F2", "rs2108622", "C", "T")
        .generate();
    assertDiplotypePairs("*2/*3", testMatchNamedAlleles(sf_definitionFile, vcfFile));
  }

  @Test
  void s2_s2s3() throws Exception {
    Path vcfFile = new TestVcfBuilder("*2/*2+*3")
        .reference("CYP4F2")
        .saveFile()
        .variation("CYP4F2", "rs3093105", "C", "C")
        .variation("CYP4F2", "rs2108622", "C", "T")
        .generate();
    assertDiplotypePairs("*2/*2 + *3", testMatchNamedAlleles(sf_definitionFile, vcfFile, true, true));
  }

  @Test
  void s2s3_s2s3() throws Exception {
    Path vcfFile = new TestVcfBuilder("*2+*3/*2+*3")
        .reference("CYP4F2")
        .saveFile()
        .variation("CYP4F2", "rs3093105", "C", "C")
        .variation("CYP4F2", "rs2108622", "T", "T")
        .generate();
    assertDiplotypePairs("*2 + *3/*2 + *3", testMatchNamedAlleles(sf_definitionFile, vcfFile, true, true));
  }
}

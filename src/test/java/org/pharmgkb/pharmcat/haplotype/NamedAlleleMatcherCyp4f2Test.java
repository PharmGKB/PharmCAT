package org.pharmgkb.pharmcat.haplotype;

import java.nio.file.Path;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.TestInfo;
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
  void s1_s4(TestInfo testInfo) throws Exception {
    Path vcfFile = new TestVcfBuilder(testInfo, "*1/*4")
        .saveFile()
        .variation("CYP4F2", "rs3093105", "C", "A")
        .variation("CYP4F2", "rs2108622", "C", "T")
        .generate();
    assertDiplotypePairs("*1/*4", testMatchNamedAlleles(sf_definitionFile, vcfFile, true));
  }

  @Test
  void s2_s4(TestInfo testInfo) throws Exception {
    Path vcfFile = new TestVcfBuilder(testInfo, "*2/*4")
        .saveFile()
        .variation("CYP4F2", "rs3093105", "C", "C")
        .variation("CYP4F2", "rs2108622", "C", "T")
        .generate();
    assertDiplotypePairs("*2/*4", testMatchNamedAlleles(sf_definitionFile, vcfFile, true, true));
  }
}

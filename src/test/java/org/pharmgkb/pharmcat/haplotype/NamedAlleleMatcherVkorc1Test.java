package org.pharmgkb.pharmcat.haplotype;

import java.nio.file.Path;
import java.util.List;
import com.google.common.collect.Lists;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.pharmgkb.common.util.PathUtils;
import org.pharmgkb.pharmcat.haplotype.model.Result;
import org.pharmgkb.pharmcat.util.DataManager;

import static org.pharmgkb.pharmcat.haplotype.NamedAlleleMatcherTest.assertDiplotypePairs;
import static org.pharmgkb.pharmcat.haplotype.NamedAlleleMatcherTest.testMatchNamedAlleles;


/**
 * Test calling for VKORC1.
 *
 * @author Lester Carter
 */
class NamedAlleleMatcherVkorc1Test {
  private Path m_definitionFile;

  @BeforeEach
  void before() {
    m_definitionFile = DataManager.DEFAULT_DEFINITION_DIR.resolve("VKORC1_translation.json");
  }

  @Test
  void vkorc1gg() throws Exception {
    // Test -1639G/-1639G

    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/VKORC1/-1639G-1639G.vcf");
    List<String> expectedMatches = Lists.newArrayList("rs9923231 reference (C)/rs9923231 reference (C)");

    Result result = testMatchNamedAlleles(m_definitionFile, vcfFile);
    assertDiplotypePairs(expectedMatches, result);
  }

  @Test
  void vkorc1ga() throws Exception {
    // Test -1639G/-1639A

    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/VKORC1/-1639G-1639A.vcf");
    List<String> expectedMatches = Lists.newArrayList("rs9923231 reference (C)/rs9923231 variant (T)");

    Result result = testMatchNamedAlleles(m_definitionFile, vcfFile);
    assertDiplotypePairs(expectedMatches, result);
  }

  @Test
  void vkorc1aa() throws Exception {
    // Test -1639A/-1639A

    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/VKORC1/-1639A-1639A.vcf");
    List<String> expectedMatches = Lists.newArrayList("rs9923231 variant (T)/rs9923231 variant (T)");

    Result result = testMatchNamedAlleles(m_definitionFile, vcfFile);
    assertDiplotypePairs(expectedMatches, result);
  }
}

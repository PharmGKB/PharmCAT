package org.pharmgkb.pharmcat.haplotype;

import java.io.IOException;
import java.nio.file.Path;
import java.util.List;
import com.google.common.collect.Lists;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.TestInfo;
import org.pharmgkb.common.util.PathUtils;
import org.pharmgkb.pharmcat.Env;
import org.pharmgkb.pharmcat.ReportableException;
import org.pharmgkb.pharmcat.TestUtils;
import org.pharmgkb.pharmcat.haplotype.model.Result;
import org.pharmgkb.pharmcat.util.DataManager;

import static org.pharmgkb.pharmcat.haplotype.NamedAlleleMatcherTest.assertDiplotypePairs;
import static org.pharmgkb.pharmcat.haplotype.NamedAlleleMatcherTest.testMatchNamedAlleles;


/**
 * Test calling for VKORC1.
 *
 * @author Lester Carter
 */
public class NamedAlleleMatcherVkorc1Test {
  private static final Path sf_definitionFile = DataManager.DEFAULT_DEFINITION_DIR.resolve("VKORC1_translation.json");
  private static Env s_env = null;

  @BeforeAll
  static void prepare() throws IOException, ReportableException {
    s_env = new Env();
    //TestUtils.setSaveTestOutput(true);
  }

  @AfterEach
  void deleteDirectory(TestInfo testInfo) {
    TestUtils.deleteTestOutputDirectory(testInfo);
  }


  @Test
  void vkorc1gg() throws Exception {
    // Test -1639G/-1639G

    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/VKORC1/-1639G-1639G.vcf");
    List<String> expectedMatches = Lists.newArrayList("rs9923231 reference (C)/rs9923231 reference (C)");

    Result result = testMatchNamedAlleles(s_env, sf_definitionFile, vcfFile);
    assertDiplotypePairs(expectedMatches, result);
  }

  @Test
  void vkorc1ga() throws Exception {
    // Test -1639G/-1639A

    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/VKORC1/-1639G-1639A.vcf");
    List<String> expectedMatches = Lists.newArrayList("rs9923231 reference (C)/rs9923231 variant (T)");

    Result result = testMatchNamedAlleles(s_env, sf_definitionFile, vcfFile);
    assertDiplotypePairs(expectedMatches, result);
  }

  @Test
  void vkorc1aa() throws Exception {
    // Test -1639A/-1639A

    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/VKORC1/-1639A-1639A.vcf");
    List<String> expectedMatches = Lists.newArrayList("rs9923231 variant (T)/rs9923231 variant (T)");

    Result result = testMatchNamedAlleles(s_env, sf_definitionFile, vcfFile);
    assertDiplotypePairs(expectedMatches, result);
  }
}

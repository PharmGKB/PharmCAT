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
 * Test calling for IFNL3.
 *
 * @author Lester Carter
 */
public class NamedAlleleMatcherIfnl3Test {
  private static final Path sf_definitionFile = DataManager.DEFAULT_DEFINITION_DIR.resolve("IFNL3_translation.json");
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
  void rs12979860CC() throws Exception {
    // Test rs12979860 CC

    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/IFNL3/rs12979860CC.vcf");
    List<String> expectedMatches = Lists.newArrayList("rs12979860 reference (C)/rs12979860 reference (C)");

    Result result = testMatchNamedAlleles(s_env, sf_definitionFile, vcfFile);
    assertDiplotypePairs(expectedMatches, result);
  }

  @Test
  void rs12979860CT() throws Exception {
    // Test rs12979860 CC

    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/IFNL3/rs12979860CT.vcf");
    List<String> expectedMatches = Lists.newArrayList("rs12979860 reference (C)/rs12979860 variant (T)");

    Result result = testMatchNamedAlleles(s_env, sf_definitionFile, vcfFile);
    assertDiplotypePairs(expectedMatches, result);
  }

  @Test
  void rs12979860TT() throws Exception {
    // Test rs12979860 CC

    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/IFNL3/rs12979860TT.vcf");
    List<String> expectedMatches = Lists.newArrayList("rs12979860 variant (T)/rs12979860 variant (T)");

    Result result = testMatchNamedAlleles(s_env, sf_definitionFile, vcfFile);
    assertDiplotypePairs(expectedMatches, result);
  }
}

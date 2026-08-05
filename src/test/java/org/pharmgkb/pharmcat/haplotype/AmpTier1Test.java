package org.pharmgkb.pharmcat.haplotype;

import java.io.IOException;
import java.nio.file.Path;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.TestInfo;
import org.pharmgkb.pharmcat.Env;
import org.pharmgkb.pharmcat.ReportableException;
import org.pharmgkb.pharmcat.TestUtils;
import org.pharmgkb.pharmcat.TestVcfBuilder;
import org.pharmgkb.pharmcat.haplotype.model.GeneCall;
import org.pharmgkb.pharmcat.haplotype.model.Result;
import org.pharmgkb.pharmcat.reporter.model.MessageAnnotation;
import org.pharmgkb.pharmcat.util.DataManager;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.pharmgkb.pharmcat.haplotype.NamedAlleleMatcherTest.testMatchNamedAlleles;


/**
 * Regression test for the "missing AMP Tier 1 position" warning ({@code DefinitionExemption.hasAmp1Positions()} /
 * {@code ResultBuilder}'s {@code missing-amp1-position} message). CYP3A4 has a single AMP1 position
 * (rs35599367 / *22) among many other, non-AMP1 positions, so omitting only that one position from the VCF leaves
 * the gene still callable while exercising the warning in isolation.
 *
 * @author Mark Woon
 */
public class AmpTier1Test {
  private static final Path sf_definitionFile = DataManager.DEFAULT_DEFINITION_DIR.resolve("CYP3A4_translation.json");
  private static Env s_env = null;

  @BeforeAll
  static void prepare() throws IOException, ReportableException {
    s_env = new Env();
  }

  @AfterEach
  void deleteDirectory(TestInfo testInfo) {
    TestUtils.deleteTestOutputDirectory(testInfo);
  }


  @Test
  void missingAmp1Position(TestInfo testInfo) throws Exception {
    Result result = testMatchNamedAlleles(s_env, sf_definitionFile, new TestVcfBuilder(testInfo)
        .reference("CYP3A4")
        .missing("CYP3A4", "rs35599367")
        .generate());

    assertEquals(1, result.getGeneCalls().size());
    GeneCall geneCall = result.getGeneCalls().get(0);
    assertEquals("CYP3A4", geneCall.getGene());

    MessageAnnotation warning = geneCall.getWarnings().stream()
        .filter(w -> "missing-amp1-position".equals(w.getName()))
        .findFirst()
        .orElse(null);
    assertTrue(warning != null, "Expected a missing-amp1-position warning");
    assertTrue(warning.getMessage().contains("chr7:99768693"),
        "Warning should cite the missing AMP Tier 1 position, but was: " + warning.getMessage());
  }
}

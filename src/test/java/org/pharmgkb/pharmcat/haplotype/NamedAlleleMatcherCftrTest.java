package org.pharmgkb.pharmcat.haplotype;

import java.nio.file.Path;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.TestInfo;
import org.pharmgkb.pharmcat.TestVcfBuilder;
import org.pharmgkb.pharmcat.util.DataManager;

import static org.pharmgkb.pharmcat.haplotype.NamedAlleleMatcherTest.assertDiplotypePairs;
import static org.pharmgkb.pharmcat.haplotype.NamedAlleleMatcherTest.testMatchNamedAlleles;


/**
 * CFTR-specific tests for {@link NamedAlleleMatcher}.
 *
 * @author Mark Woon
 */
class NamedAlleleMatcherCftrTest {
  private static final Path sf_definitionFile = DataManager.DEFAULT_DEFINITION_DIR.resolve("CFTR_translation.json");


  @Test
  void reference_Reference(TestInfo testInfo) throws Exception {
    assertDiplotypePairs("ivacaftor non-responsive CFTR sequence/ivacaftor non-responsive CFTR sequence",
        testMatchNamedAlleles(sf_definitionFile, new TestVcfBuilder(testInfo, "ref/ref")
            .reference("CFTR")
            .generate()));
  }


  @Test
  void G1244E_reference(TestInfo testInfo) throws Exception {
    assertDiplotypePairs("G1244E/ivacaftor non-responsive CFTR sequence",
        testMatchNamedAlleles(sf_definitionFile, new TestVcfBuilder(testInfo, "G1244E/ref")
            .variation("CFTR", "rs267606723", "G", "A")
            .generate()));
  }
}

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
 * JUnit test for {@link NamedAlleleMatcher#callDiplotypes(MatchData, boolean)}.
 *
 * @author Lester Carter
 */
class NamedAlleleMatcherDpydTest {
  private Path m_definitionFile;

  @BeforeEach
  void before() {
    m_definitionFile = DataManager.DEFAULT_DEFINITION_DIR.resolve("DPYD_translation.json");
  }


  @Test
  void dpyds1s1() throws Exception {
    // Test *1/*1 - contains a del

    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/DPYD/novariant.vcf");
    List<String> expectedMatches = Lists.newArrayList("Reference/Reference");

    Result result = testMatchNamedAlleles(m_definitionFile, vcfFile);
    assertDiplotypePairs(expectedMatches, result);
  }


  @Test
  void dpyds2aRs67376798A() throws Exception {
    // Test *2a/Rs67376798A

    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/DPYD/s2aRs67376798A.vcf");
    List<String> expectedMatches = Lists.newArrayList("c.1905+1G>A (*2A)/c.2846A>T");

    Result result = testMatchNamedAlleles(m_definitionFile, vcfFile);
    assertDiplotypePairs(expectedMatches, result);
  }

  @Test
  void dpyds1s2b() throws Exception {
    // Test *1/*2b - however can't be distinguished from *2A/*5

    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/DPYD/s1s2b.vcf");
    List<String> expectedMatches = Lists.newArrayList("c.1627A>G (*5)/c.1905+1G>A (*2A)");

    Result result = testMatchNamedAlleles(m_definitionFile, vcfFile);
    assertDiplotypePairs(expectedMatches, result);
  }
}

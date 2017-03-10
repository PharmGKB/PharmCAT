package org.pharmgkb.pharmcat.haplotype;

import java.nio.file.Path;
import java.util.List;
import com.google.common.collect.Lists;
import org.junit.Before;
import org.junit.Test;
import org.pharmgkb.common.util.PathUtils;
import org.pharmgkb.pharmcat.definition.DefinitionManager;
import org.pharmgkb.pharmcat.haplotype.model.Result;

import static org.pharmgkb.pharmcat.haplotype.NamedAlleleMatcherTest.assertDiplotypePairs;
import static org.pharmgkb.pharmcat.haplotype.NamedAlleleMatcherTest.testMatchNamedAlleles;


/**
 * JUnit test for {@link NamedAlleleMatcher#callDiplotypes(MatchData, boolean)}.
 *
 * @author Lester Carter
 */
public class NamedAlleleMatcherIfnl3Test {
  private Path m_definitionFile;

  @Before
  public void before() throws Exception {
    m_definitionFile = DefinitionManager.DEFAULT_DEFINITION_DIR.resolve("IFNL3_translation.json");
  }

  @Test
  public void rs12979860CC() throws Exception {
    // Test rs12979860 CC

    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/IFNL3/rs12979860CC.vcf");
    List<String> expectedMatches = Lists.newArrayList("rs12979860C/rs12979860C");

    Result result = testMatchNamedAlleles(m_definitionFile, vcfFile);
    assertDiplotypePairs(expectedMatches, result);
  }

  @Test
  public void rs12979860CT() throws Exception {
    // Test rs12979860 CC

    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/IFNL3/rs12979860CT.vcf");
    List<String> expectedMatches = Lists.newArrayList("rs12979860C/rs12979860T");

    Result result = testMatchNamedAlleles(m_definitionFile, vcfFile);
    assertDiplotypePairs(expectedMatches, result);
  }

  @Test
  public void rs12979860TT() throws Exception {
    // Test rs12979860 CC

    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/IFNL3/rs12979860TT.vcf");
    List<String> expectedMatches = Lists.newArrayList("rs12979860T/rs12979860T");

    Result result = testMatchNamedAlleles(m_definitionFile, vcfFile);
    assertDiplotypePairs(expectedMatches, result);
  }
}

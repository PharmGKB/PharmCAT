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
public class NamedAlleleMatcherCyp3a5Test {
  private Path m_definitionFile;

  @Before
  public void before() throws Exception {
    m_definitionFile = DefinitionManager.DEFAULT_DEFINITION_DIR.resolve("CYP3A5_translation.json");
  }

  @Test
  public void cyp3a5s3s9() throws Exception {
    // Test *3/*9.

    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/cyp3a5/s3s9.vcf");
    List<String> expectedMatches = Lists.newArrayList("*3/*9");

    Result result = testMatchNamedAlleles(m_definitionFile, vcfFile);
    assertDiplotypePairs(expectedMatches, result);
  }

  @Test
  public void cyp3a5s1s7() throws Exception {
    // Test *1/*7.  Het in rs41303343 position, a single insertion.

    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/cyp3a5/s1s7.vcf");
    List<String> expectedMatches = Lists.newArrayList("*1/*7");

    Result result = testMatchNamedAlleles(m_definitionFile, vcfFile);
    assertDiplotypePairs(expectedMatches, result);
  }


  @Test
  public void cyp3a5s3s9Homozygous() throws Exception {
    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/cyp3a5/s3s9-homozygous.vcf");
    List<String> expectedMatches = Lists.newArrayList("*3/*9");

    Result result = testMatchNamedAlleles(m_definitionFile, vcfFile);
    assertDiplotypePairs(expectedMatches, result);
  }
}

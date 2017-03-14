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
public class NamedAlleleMatcherCftrTest {
  private Path m_definitionFile;

  @Before
  public void before() throws Exception {
    m_definitionFile = DefinitionManager.DEFAULT_DEFINITION_DIR.resolve("CFTR_translation.json");
  }


  @Test
  public void cftrReferenceReference() throws Exception {
    // Test reference
    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/cftr/refref.vcf");
    List<String> expectedMatches = Lists.newArrayList("Reference/Reference");

    Result result = testMatchNamedAlleles(m_definitionFile, vcfFile);
    assertDiplotypePairs(expectedMatches, result);
  }


  @Test
  public void cftrF508delF508del() throws Exception {
    // Test F508del(TCT)/F508del(TCT)
    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/cftr/F508delF508del.vcf");
    List<String> expectedMatches = Lists.newArrayList("F508del(TCT)/F508del(TCT)");

    Result result = testMatchNamedAlleles(m_definitionFile, vcfFile);
    assertDiplotypePairs(expectedMatches, result);
  }

}

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
class NamedAlleleMatcherCftrTest {
  private Path m_definitionFile;

  @BeforeEach
  void before() {
    m_definitionFile = DataManager.DEFAULT_DEFINITION_DIR.resolve("CFTR_translation.json");
  }


  @Test
  void cftrReferenceReference() throws Exception {
    // Test reference
    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/cftr/refref.vcf");
    List<String> expectedMatches = Lists.newArrayList("ivacaftor non-responsive CFTR sequence/ivacaftor non-responsive CFTR sequence");

    Result result = testMatchNamedAlleles(m_definitionFile, vcfFile);
    assertDiplotypePairs(expectedMatches, result);
  }


  @Test
  void cftrF508delF508del() throws Exception {
    // Test F508del(TCT)/F508del(TCT)
    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/cftr/F508delF508del.vcf");
    List<String> expectedMatches = Lists.newArrayList("ivacaftor non-responsive CFTR sequence/ivacaftor non-responsive CFTR sequence");

    Result result = testMatchNamedAlleles(m_definitionFile, vcfFile);
    assertDiplotypePairs(expectedMatches, result);
  }


  @Test
  void G1244Eref() throws Exception {
    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/cftr/G1244Eref.vcf");
    List<String> expectedMatches = Lists.newArrayList("G1244E/ivacaftor non-responsive CFTR sequence");

    Result result = testMatchNamedAlleles(m_definitionFile, vcfFile);
    assertDiplotypePairs(expectedMatches, result);
  }


  @Test
  void G1244EF508del() throws Exception {
    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/cftr/G1244EF508del.vcf");
    List<String> expectedMatches = Lists.newArrayList("G1244E/ivacaftor non-responsive CFTR sequence");

    Result result = testMatchNamedAlleles(m_definitionFile, vcfFile);
    assertDiplotypePairs(expectedMatches, result);
  }


  @Test
  void G551DG542X() throws Exception {
    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/cftr/G551DG542X.vcf");
    List<String> expectedMatches = Lists.newArrayList("G551D/ivacaftor non-responsive CFTR sequence");

    Result result = testMatchNamedAlleles(m_definitionFile, vcfFile);
    assertDiplotypePairs(expectedMatches, result);
  }


}

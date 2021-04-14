package org.pharmgkb.pharmcat.haplotype;

import java.nio.file.Path;
import java.util.List;
import java.util.SortedSet;
import com.google.common.collect.Lists;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.pharmgkb.common.util.PathUtils;
import org.pharmgkb.pharmcat.haplotype.model.GeneCall;
import org.pharmgkb.pharmcat.haplotype.model.Result;
import org.pharmgkb.pharmcat.haplotype.model.Variant;
import org.pharmgkb.pharmcat.util.DataManager;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.pharmgkb.pharmcat.haplotype.NamedAlleleMatcherTest.assertDiplotypePairs;
import static org.pharmgkb.pharmcat.haplotype.NamedAlleleMatcherTest.testMatchNamedAlleles;


/**
 * JUnit test for {@link NamedAlleleMatcher#callDiplotypes(MatchData, boolean)}.
 *
 * @author Lester Carter
 */
class NamedAlleleMatcherCyp2c9Test {
  private Path m_definitionFile;

  @BeforeEach
  void before() {
    m_definitionFile = DataManager.DEFAULT_DEFINITION_DIR.resolve("CYP2C9_translation.json");
  }

  @Test
  void cyp2c9s1s1() throws Exception {
    // Test *1/*1

    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/cyp2c9/s1s1.vcf");
    List<String> expectedMatches = Lists.newArrayList("*1/*1");

    Result result = testMatchNamedAlleles(m_definitionFile, vcfFile);
    assertDiplotypePairs(expectedMatches, result);
  }

  @Test
  void cyp2c9s2s3() throws Exception {
    // Test *2/*3

    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/cyp2c9/s2s3.vcf");
    List<String> expectedMatches = Lists.newArrayList("*2/*3");

    Result result = testMatchNamedAlleles(m_definitionFile, vcfFile);
    assertDiplotypePairs(expectedMatches, result);
  }

  @Test
  void cyp2c9s2s24() throws Exception {

    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/cyp2c9/s2s24.vcf");
    List<String> expectedMatches = Lists.newArrayList("*2/*24");

    Result result = testMatchNamedAlleles(m_definitionFile, vcfFile);
    assertDiplotypePairs(expectedMatches, result);
  }

  @Test
  void cyp2c9s2s24Only() throws Exception {

    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/cyp2c9/s2s24only.vcf");
    List<String> expectedMatches = Lists.newArrayList(); // no expected match

    Result result = testMatchNamedAlleles(m_definitionFile, vcfFile);
    assertDiplotypePairs(expectedMatches, result);
  }

  @Test
  void cyp2c9s24s24() throws Exception {
    // Test *24/*24 with hom at rs749060448

    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/cyp2c9/s24s24.vcf");
    List<String> expectedMatches = Lists.newArrayList("*24/*24");

    Result result = testMatchNamedAlleles(m_definitionFile, vcfFile);
    assertDiplotypePairs(expectedMatches, result);
  }

  @Test
  void cyp2c9s2s25() throws Exception {
    // Test *2/*25

    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/cyp2c9/s2s25.vcf");
    List<String> expectedMatches = Lists.newArrayList("*2/*25");

    Result result = testMatchNamedAlleles(m_definitionFile, vcfFile);
    assertDiplotypePairs(expectedMatches, result);
  }

  /**
   * This tests what happens when a homozygous allele is mixed with a heterozygous allele. It should result in a
   * non-call since the matcher cannot match to a single definition.
   *
   * This also asserts that no missing or mismatched data exists to ensure the mismatch is not due to invalid input.
   */
  @Test
  void cyp2c9s24s2s24() throws Exception {
    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/cyp2c9/s24s2pluss24.vcf");

    Result result = testMatchNamedAlleles(m_definitionFile, vcfFile, true, false, true, true);
    GeneCall call = result.getGeneCalls().stream()
        .filter(c -> c.getGene().equals("CYP2C9")).findFirst()
        .orElseThrow(() -> new RuntimeException("No gene call found"));
    assertEquals(0, call.getDiplotypes().size());
    assertEquals(0, call.getMatchData().getMissingPositions().size());
    assertEquals(0, call.getMatchData().getMismatchedPositions().size());
  }


  @Test
  void testExtraPosition() throws Exception {

    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/cyp2c9/s1s1.vcf");

    DefinitionReader definitionReader = new DefinitionReader();
    definitionReader.read(m_definitionFile);
    definitionReader.readExemptions(DataManager.DEFAULT_DEFINITION_DIR.resolve(DataManager.EXEMPTIONS_JSON_FILE_NAME));

    NamedAlleleMatcher namedAlleleMatcher = new NamedAlleleMatcher(definitionReader, true, false);
    Result result = namedAlleleMatcher.call(vcfFile);
    SortedSet<Variant> extraPositions = result.getGeneCalls().get(0).getVariantsOfInterest();
    assertEquals(1, extraPositions.size());
    assertEquals("rs12777823", extraPositions.first().getRsid());
    assertEquals("G/A", extraPositions.first().getVcfCall());
  }
}

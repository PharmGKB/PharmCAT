package org.pharmgkb.pharmcat.haplotype;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.SortedSet;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.TestInfo;
import org.pharmgkb.pharmcat.Env;
import org.pharmgkb.pharmcat.TestVcfBuilder;
import org.pharmgkb.pharmcat.VcfFile;
import org.pharmgkb.pharmcat.definition.DefinitionReader;
import org.pharmgkb.pharmcat.haplotype.model.GeneCall;
import org.pharmgkb.pharmcat.haplotype.model.Result;
import org.pharmgkb.pharmcat.haplotype.model.Variant;
import org.pharmgkb.pharmcat.util.DataManager;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.pharmgkb.pharmcat.haplotype.NamedAlleleMatcherTest.assertDiplotypePairs;
import static org.pharmgkb.pharmcat.haplotype.NamedAlleleMatcherTest.testMatchNamedAlleles;


/**
 * CYP2C9-specific tests for {@link NamedAlleleMatcher}.
 *
 * @author Mark Woon
 */
class NamedAlleleMatcherCyp2c9Test {
  private final Path sf_definitionFile = DataManager.DEFAULT_DEFINITION_DIR.resolve("CYP2C9_translation.json");


  @Test
  void s1_s1(TestInfo testInfo) throws Exception {
    assertDiplotypePairs("*1/*1", testMatchNamedAlleles(sf_definitionFile,
        new TestVcfBuilder(testInfo, "*1/*1")
            .reference("CYP2C9")
            .generate()));
  }

  @Test
  void s2_s3(TestInfo testInfo) throws Exception {
    assertDiplotypePairs("*2/*3", testMatchNamedAlleles(sf_definitionFile,
        new TestVcfBuilder(testInfo, "*2/*3")
            .variation("CYP2C9", "rs1799853", "C", "T")
            .variation("CYP2C9", "rs1057910", "A", "C")
            .generate()));
  }

  @Test
  void s2_s24(TestInfo testInfo) throws Exception {
    assertDiplotypePairs("*2/*24", testMatchNamedAlleles(sf_definitionFile,
        new TestVcfBuilder(testInfo, "*2/*24")
            .variation("CYP2C9", "rs1799853", "C", "T")
            .variation("CYP2C9", "rs749060448", "A", "G")
            .generate()));
  }

  @Test
  void s2_s24Only(TestInfo testInfo) throws Exception {
    // no expected match
    assertDiplotypePairs(new ArrayList<>(), testMatchNamedAlleles(sf_definitionFile,
        new TestVcfBuilder(testInfo, "*2/*24 only")
        .variation("CYP2C9", "rs1799853", "C", "T")
        .variation("CYP2C9", "rs749060448", "A", "A")
        .generate()));
  }

  @Test
  void s24_s24(TestInfo testInfo) throws Exception {
    assertDiplotypePairs("*24/*24", testMatchNamedAlleles(sf_definitionFile,
        new TestVcfBuilder(testInfo, "*24/*24")
            .variation("CYP2C9", "rs749060448", "A", "A")
            .generate()));
  }

  @Test
  void s2_s25(TestInfo testInfo) throws Exception {
    assertDiplotypePairs("*2/*25", testMatchNamedAlleles(sf_definitionFile,
        new TestVcfBuilder(testInfo, "*2/*25")
            .variation("CYP2C9", "rs1799853", "C", "T")
            .variation("CYP2C9", "rs1304490498", "AGAAATGGAA", "delAGAAATGGAA")
            .generate()));
  }

  /**
   * This tests what happens when a homozygous allele is mixed with a heterozygous allele. It should result in a
   * non-call since the matcher cannot match to a single definition.
   * <p>
   * This also asserts that no missing or mismatched data exists to ensure the mismatch is not due to invalid input.
   */
  @Test
  void s24_s2s24(TestInfo testInfo) throws Exception {
    Path vcfFile = new TestVcfBuilder(testInfo, "*24/*2+*24")
        .variation("CYP2C9", "rs1799853", "C", "T")
        .variation("CYP2C9", "rs749060448", "A", "A")
        .generate();

    Result result = testMatchNamedAlleles(sf_definitionFile, vcfFile, false);
    List<GeneCall> calls = result.getGeneCalls().stream()
        .filter(c -> c.getGene().equals("CYP2C9"))
        .toList();
    assertEquals(1, calls.size(), "No gene call found");
    GeneCall call = calls.get(0);
    System.out.println(call.getDiplotypes());
    assertEquals(0, call.getDiplotypes().size());
    assertEquals(0, call.getMatchData().getMissingPositions().size());
    assertEquals(0, call.getMatchData().getPositionsWithUndocumentedVariations().size());
  }


  @Test
  void extraPosition(TestInfo testInfo) throws Exception {
    Path vcfFile = new TestVcfBuilder(testInfo, "*1/*1")
        .reference("CYP2C9")
        .extraPosition("CYP2C9", "rs12777823", "G", "A")
        .generate();

    DefinitionReader definitionReader = new DefinitionReader(sf_definitionFile, DataManager.DEFAULT_EXEMPTIONS_FILE);

    NamedAlleleMatcher namedAlleleMatcher = new NamedAlleleMatcher(new Env(), definitionReader, false, false, true);
    Result result = namedAlleleMatcher.call(new VcfFile(vcfFile, false), null);
    SortedSet<Variant> extraPositions = result.getGeneCalls().get(0).getVariantsOfInterest();
    assertEquals(1, extraPositions.size());
    assertEquals("rs12777823", extraPositions.first().getRsid());
    assertEquals("G/A", extraPositions.first().getVcfCall());
  }
}

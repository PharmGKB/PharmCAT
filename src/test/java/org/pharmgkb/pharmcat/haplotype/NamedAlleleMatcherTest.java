package org.pharmgkb.pharmcat.haplotype;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Objects;
import java.util.Optional;
import java.util.Set;
import java.util.TreeSet;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import com.google.common.base.Joiner;
import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.TestInfo;
import org.pharmgkb.common.util.PathUtils;
import org.pharmgkb.pharmcat.DiplotypeUtils;
import org.pharmgkb.pharmcat.Env;
import org.pharmgkb.pharmcat.TestVcfBuilder;
import org.pharmgkb.pharmcat.VcfFile;
import org.pharmgkb.pharmcat.definition.DefinitionReader;
import org.pharmgkb.pharmcat.definition.model.NamedAllele;
import org.pharmgkb.pharmcat.haplotype.model.BaseMatch;
import org.pharmgkb.pharmcat.haplotype.model.CombinationMatch;
import org.pharmgkb.pharmcat.haplotype.model.DiplotypeMatch;
import org.pharmgkb.pharmcat.haplotype.model.GeneCall;
import org.pharmgkb.pharmcat.haplotype.model.Result;
import org.pharmgkb.pharmcat.util.DataManager;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.contains;
import static org.hamcrest.core.IsEqual.equalTo;
import static org.junit.jupiter.api.Assertions.*;


/**
 * JUnit test for {@link NamedAlleleMatcher}.
 *
 * @author Mark Woon
 */
class NamedAlleleMatcherTest {


  /**
   * Helper method for running {@link NamedAlleleMatcher} using:
   * <ul>
   *   <li>{@code topCandidateOnly} = false</li>
   *   <li>{@code showUnmatched} = true</li>
   *   <li>{@code withExemptions} = true</li>
   * </ul>
   *
   */
  static Result testMatchNamedAlleles(Path tsvFile, Path vcfFile) throws Exception {
    return testMatchNamedAlleles(tsvFile, vcfFile, false, false, true, true);
  }

  static Result testMatchNamedAlleles(Path tsvFile, Path vcfFile, boolean topCandidateOnly)
      throws Exception {
    return testMatchNamedAlleles(tsvFile, vcfFile, false, topCandidateOnly, true, true);
  }

  static Result testMatchNamedAlleles(Path tsvFile, Path vcfFile, boolean findCombination, boolean topCandidateOnly)
      throws Exception {
    return testMatchNamedAlleles(tsvFile, vcfFile, findCombination, topCandidateOnly, true, true);
  }

  /**
   * Helper method for running {@link NamedAlleleMatcher}.
   * This is used by the more specific gene tests.
   */
  static Result testMatchNamedAlleles(Path definitionFile, Path vcfFile, boolean topCandidateOnly,
      boolean showUnmatched, boolean withExemptions) throws Exception {
    return testMatchNamedAlleles(definitionFile, vcfFile, false, topCandidateOnly, showUnmatched,
        withExemptions);
  }

  static Result testMatchNamedAlleles(Path definitionFile, Path vcfFile, boolean findCombinations,
      boolean topCandidateOnly, boolean showUnmatched, boolean withExemptions) throws Exception {

    DefinitionReader definitionReader = new DefinitionReader(definitionFile,
        withExemptions ? DataManager.DEFAULT_EXEMPTIONS_FILE : null);

    NamedAlleleMatcher namedAlleleMatcher = new NamedAlleleMatcher(new Env(), definitionReader, findCombinations,
        topCandidateOnly, true);
    Result result = namedAlleleMatcher.call(new VcfFile(vcfFile), null);

    // print
    new ResultSerializer()
        .alwaysShowUnmatchedHaplotypes(showUnmatched)
        .toHtml(result, vcfFile.getParent().resolve(PathUtils.getBaseFilename(vcfFile) + ".html"))
        .toJson(result, vcfFile.getParent().resolve(PathUtils.getBaseFilename(vcfFile) + ".json"));

    return result;
  }


  /**
   * Checks that the list of diplotype matches are what we expect.
   *
   * @param expectedPair the expected diplotype in "*1/*2" format
   * @param result the {@link NamedAlleleMatcher} results
 0  */
  static void assertDiplotypePairs(String expectedPair, Result result) {
    assertDiplotypePairs(Lists.newArrayList(expectedPair), result);
  }

  /**
   * Checks that the list of diplotype matches are what we expect.
   *
   * @param expectedPairs the set of expected diplotypes in "*1/*2" format
   * @param result the {@link NamedAlleleMatcher} results
   */
  static void assertDiplotypePairs(List<String> expectedPairs, Result result) {

    Preconditions.checkNotNull(expectedPairs);
    Preconditions.checkNotNull(result);

    List<String> pairs = new ArrayList<>();
    StringBuilder builder = new StringBuilder();
    if (result.getGeneCalls().size() > 0) {
      for (DiplotypeMatch dm : result.getGeneCalls().get(0).getDiplotypes()) {
        pairs.add(dm.getName());
        if (builder.length() > 0) {
          builder.append(", ");
        }
        builder.append(dm.getName())
            .append(" (")
            .append(dm.getScore())
            .append(")");
      }
    }

    if (expectedPairs.size() != pairs.size() || !expectedPairs.equals(pairs)) {
      System.out.println("Expected: [" + Joiner.on(", ").join(expectedPairs) + "]");
      System.out.println("Got:      " + pairs);
      System.out.println("Scores:   " + builder);
      fail("Did not get expected matches");
    }
  }

  static void assertPhasedOutput(String expected, Result result) {
    List<String> pairs = new ArrayList<>();
    GeneCall geneCall = result.getGeneCalls().get(0);

    assertTrue(geneCall.isPhased(), "Provided sample is not phased so this test is improper");

    for (DiplotypeMatch dm : result.getGeneCalls().get(0).getDiplotypes()) {
      pairs.add(dm.getName());
    }

    String testResult = pairs.stream().reduce(DiplotypeUtils.PhasedReducer).orElseThrow(RuntimeException::new);

    assertEquals(expected, testResult, "Phased output wasn't expected");
  }


  @Test
  void testCall() throws Exception {

    Path vcfFile  = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/haplotyper.vcf");
    Path jsonFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/haplotyper.json");

    DefinitionReader definitionReader = new DefinitionReader(jsonFile, null);

    NamedAlleleMatcher namedAlleleMatcher = new NamedAlleleMatcher(new Env(), definitionReader, true, true, true);
    Result result = namedAlleleMatcher.call(new VcfFile(vcfFile), null);
    Set<DiplotypeMatch> pairs = result.getGeneCalls().get(0).getDiplotypes();
    assertNotNull(pairs);
    assertEquals(1, pairs.size());
    assertEquals("*1/*2", pairs.iterator().next().getName());
  }


  /**
   * This breaks down the main code path that {@link #testCall()} runs to simplify testing smaller chunks at a time.
   */
  @Test
  void testCallDiplotypePath() throws Exception {

    Path vcfPath  = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/haplotyper.vcf");
    Path jsonFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/haplotyper.json");
    String gene = "CYP3A5";
    String sampleId = "NA12878";

    DefinitionReader definitionReader = new DefinitionReader(jsonFile, null);

    VcfFile vcfFile = new VcfFile(vcfPath, false);
    VcfReader vcfReader = vcfFile.getReader(definitionReader, sampleId, false);

    // grab SampleAlleles for all positions related to current gene
    MatchData data = new MatchData(sampleId, gene, vcfReader.getAlleleMap(), definitionReader.getPositions(gene), null, null);
    assertEquals(3, data.getNumSampleAlleles());
    assertEquals(0, data.getMissingPositions().size());
    // handle missing positions of interest in sample
    data.marshallHaplotypes("TEST", definitionReader.getHaplotypes(gene), false);
    assertEquals(3, data.getPositions().length);
    assertEquals(2, data.getHaplotypes().size());

    // get all permutations of sample at positions of interest
    Set<String> permutations = Sets.newHashSet(
        "1:C;2:C;4:TG;",
        "1:T;2:CT;4:T;"
    );
    data.generateSamplePermutations();
    assertThat(data.getPermutations(), equalTo(permutations));

    List<DiplotypeMatch> pairs = new DiplotypeMatcher(data).compute(false);
    assertNotNull(pairs);
    assertEquals(1, pairs.size());
    assertEquals("*1/*2", pairs.get(0).getName());
  }


  @Test
  void testMismatchedRefAlleleWarnings() throws Exception {

    Path definitionFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/NamedAlleleMatcher-mismatchedRefAllele.json");
    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/NamedAlleleMatcher-mismatchedRefAllele.vcf");

    DefinitionReader definitionReader = new DefinitionReader(definitionFile, null);
    NamedAlleleMatcher namedAlleleMatcher = new NamedAlleleMatcher(new Env(), definitionReader, true, true, true);
    Result result = namedAlleleMatcher.call(new VcfFile(vcfFile), null);
    assertNotNull(result.getVcfWarnings());

    for (String key : result.getVcfWarnings().keySet()) {
      System.out.println(key);
      System.out.println("\t" + result.getVcfWarnings().get(key));
    }

    assertEquals(2, result.getVcfWarnings().size());

    assertNotNull(          result.getVcfWarnings().get("chr10:94942205"));
    assertEquals(1, result.getVcfWarnings().get("chr10:94942205").size());
    assertTrue(result.getVcfWarnings().get("chr10:94942205").iterator().next()
        .contains("does not match expected reference"));

    assertNotNull(          result.getVcfWarnings().get("chr10:94949281"));
    assertEquals(1, result.getVcfWarnings().get("chr10:94949281").size());
    assertTrue(result.getVcfWarnings().get("chr10:94949281").iterator().next()
        .contains("does not match expected reference"));
  }


  @Test
  void testCyp2d6() throws Exception {

    Path definitionFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/NamedAlleleMatcher-cyp2d6.json");
    Path vcfPath = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/NamedAlleleMatcher-cyp2d6.vcf");

    DefinitionReader definitionReader = new DefinitionReader(definitionFile, null);
    VcfFile vcfFile = new VcfFile(vcfPath);

    NamedAlleleMatcher naNoCyp2d6 = new NamedAlleleMatcher(new Env(), definitionReader, false, false, false);
    Result result = naNoCyp2d6.call(vcfFile, null);
    assertEquals(0, result.getVcfWarnings().size());
    assertEquals(0, result.getGeneCalls().size());

    NamedAlleleMatcher naWithCyp2d6 = new NamedAlleleMatcher(new Env(), definitionReader, false, false, true);
    result = naWithCyp2d6.call(vcfFile, null);
    assertEquals(0, result.getVcfWarnings().size());
    assertEquals(1, result.getGeneCalls().size());
  }


  @Test
  void testWobbleScoring(TestInfo testInfo) throws Exception {
    Path definitionFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/NamedAlleleMatcher-cyp2d6.json");
    // without wobble scoring, if wobble is reference, this will only return *1/*2
    Path vcfFile = new TestVcfBuilder(testInfo, "*1/*12")
        .withDefinition(definitionFile)
        .variation("CYP2D6", "rs1135840", "C", "G")
        .variation("CYP2D6", "rs16947", "G", "A")
        // wobble on *12
        .variation("CYP2D6", "rs28371710", "C", "C")
        .variation("CYP2D6", "rs1058164", "G", "G")
        .missing("CYP2D6", "rs5030862")
        .phased()
        .generate();

    DefinitionReader definitionReader = new DefinitionReader(definitionFile, null);
    NamedAlleleMatcher namedAlleleMatcher = new NamedAlleleMatcher(new Env(), definitionReader, false, true, true);
    Result result = namedAlleleMatcher.call(new VcfFile(vcfFile), null);
    assertEquals(0, result.getVcfWarnings().size());
    assertEquals(1, result.getGeneCalls().size());

    GeneCall geneCall = result.getGeneCalls().get(0);
    List<String> matches = printMatches(geneCall);
    assertEquals(2, matches.size());
    assertThat(matches, contains("*1/*2", "*1/*12"));
  }


  @Test
  void testCombinationBaseline() throws Exception {
    Path definitionFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/NamedAlleleMatcher-combination.json");
    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/NamedAlleleMatcher-combinationBaseline.vcf");

    DefinitionReader definitionReader = new DefinitionReader(definitionFile, null);
    NamedAlleleMatcher namedAlleleMatcher = new NamedAlleleMatcher(new Env(), definitionReader, true, true, false);
    Result result = namedAlleleMatcher.call(new VcfFile(vcfFile), null);
    assertEquals(0, result.getVcfWarnings().size());
    assertEquals(1, result.getGeneCalls().size());

    GeneCall geneCall = result.getGeneCalls().get(0);
    assertEquals(1, geneCall.getDiplotypes().size());
    assertEquals("*1/*1", geneCall.getDiplotypes().iterator().next().getName());
  }

  @Test
  void testCombinationHomozygous(TestInfo testInfo) throws Exception {
    Path definitionFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/NamedAlleleMatcher-cyp2b6.json");
    Path vcfFile = new TestVcfBuilder(testInfo, "*2 + *5/*2 + *5")
        .withDefinition(definitionFile)
        // *2
        .variation("CYP2B6", "rs8192709", "T", "T")
        // *5
        .variation("CYP2B6", "rs3211371", "T", "T")
        .generate();

    DefinitionReader definitionReader = new DefinitionReader(definitionFile, null);
    NamedAlleleMatcher namedAlleleMatcher = new NamedAlleleMatcher(new Env(), definitionReader, true, true, false);
    Result result = namedAlleleMatcher.call(new VcfFile(vcfFile), null);
    assertEquals(0, result.getVcfWarnings().size());
    assertEquals(1, result.getGeneCalls().size());

    GeneCall geneCall = result.getGeneCalls().get(0);
    assertEquals(1, geneCall.getDiplotypes().size());
    assertEquals("[*2 + *5]/[*2 + *5]", geneCall.getDiplotypes().iterator().next().getName());
  }

  @Test
  void testCombinationPhased() throws Exception {
    Path definitionFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/NamedAlleleMatcher-combination.json");
    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/NamedAlleleMatcher-combinationPhased.vcf");

    DefinitionReader definitionReader = new DefinitionReader(definitionFile, null);
    NamedAlleleMatcher namedAlleleMatcher = new NamedAlleleMatcher(new Env(), definitionReader, true, true, false);
    Result result = namedAlleleMatcher.call(new VcfFile(vcfFile), null);
    assertEquals(0, result.getVcfWarnings().size());
    assertEquals(1, result.getGeneCalls().size());

    GeneCall geneCall = result.getGeneCalls().get(0);
    printMatches(geneCall);
    assertEquals(1, geneCall.getDiplotypes().size());
    DiplotypeMatch dm = geneCall.getDiplotypes().iterator().next();
    assertEquals("*1/[*6 + *27 + *28 + *80]", dm.getName());
    assertFalse(dm.getHaplotype1().getHaplotype().isCombination());
    assertNotNull(dm.getHaplotype2());
    assertTrue(dm.getHaplotype2().getHaplotype().isCombination());
  }

  @Test
  void testCombinationUnphased() throws Exception {
    Path definitionFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/NamedAlleleMatcher-combination.json");
    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/NamedAlleleMatcher-combinationUnphased.vcf");

    DefinitionReader definitionReader = new DefinitionReader(definitionFile, null);
    NamedAlleleMatcher namedAlleleMatcher = new NamedAlleleMatcher(new Env(), definitionReader, true, true, false);
    Result result = namedAlleleMatcher.call(new VcfFile(vcfFile), null);
    assertEquals(0, result.getVcfWarnings().size());
    assertEquals(1, result.getGeneCalls().size());

    GeneCall geneCall = result.getGeneCalls().get(0);
    printMatches(geneCall);
    assertEquals(1, geneCall.getDiplotypes().size());
    DiplotypeMatch dm = geneCall.getDiplotypes().iterator().next();
    assertEquals("*1/[*6 + *27 + *28 + *80]", dm.getName());
    assertFalse(dm.getHaplotype1().getHaplotype().isCombination());
    assertNotNull(dm.getHaplotype2());
    assertTrue(dm.getHaplotype2().getHaplotype().isCombination());
  }

  @Test
  void testPartialWithCombination() throws Exception {
    Path definitionFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/NamedAlleleMatcher-combination.json");
    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/NamedAlleleMatcher-partialWithCombination.vcf");

    DefinitionReader definitionReader = new DefinitionReader(definitionFile, null);
    NamedAlleleMatcher namedAlleleMatcher = new NamedAlleleMatcher(new Env(), definitionReader, true, true, false);
    Result result = namedAlleleMatcher.call(new VcfFile(vcfFile), null);
    assertEquals(1, result.getVcfWarnings().size());
    assertEquals(1, result.getGeneCalls().size());

    GeneCall geneCall = result.getGeneCalls().get(0);
    printMatches(geneCall);
    assertEquals(1, geneCall.getDiplotypes().size());
    DiplotypeMatch dm = geneCall.getDiplotypes().iterator().next();
    assertEquals("*1/[*6 + *28 + g.233760973C>T]", dm.getName());
    assertFalse(dm.getHaplotype1().getHaplotype().isCombination());
    assertFalse(dm.getHaplotype1().getHaplotype().isPartial());
    assertNotNull(dm.getHaplotype2());
    assertTrue(dm.getHaplotype2().getHaplotype().isCombination());
    assertTrue(dm.getHaplotype2().getHaplotype().isPartial());
  }

  @Test
  void testPartial() throws Exception {
    Path definitionFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/NamedAlleleMatcher-combination.json");
    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/NamedAlleleMatcher-partial.vcf");

    DefinitionReader definitionReader = new DefinitionReader(definitionFile, null);
    NamedAlleleMatcher namedAlleleMatcher = new NamedAlleleMatcher(new Env(), definitionReader, true, true, false);
    Result result = namedAlleleMatcher.call(new VcfFile(vcfFile), null);
    assertEquals(1, result.getVcfWarnings().size());
    assertEquals(1, result.getGeneCalls().size());

    GeneCall geneCall = result.getGeneCalls().get(0);
    printMatches(geneCall);
    assertEquals(1, geneCall.getDiplotypes().size());
    DiplotypeMatch dm = geneCall.getDiplotypes().iterator().next();
    assertEquals("*6/[*6 + g.233760973C>T]", dm.getName());
    assertFalse(dm.getHaplotype1().getHaplotype().isCombination());
    assertFalse(dm.getHaplotype1().getHaplotype().isPartial());
    assertNotNull(dm.getHaplotype2());
    assertFalse(dm.getHaplotype2().getHaplotype().isCombination());
    assertTrue(dm.getHaplotype2().getHaplotype().isPartial());
  }

  @Test
  void testPartial2Phased() throws Exception {
    Path definitionFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/NamedAlleleMatcher-cyp2c19.json");
    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/NamedAlleleMatcher-partial2Phased.vcf");

    DefinitionReader definitionReader = new DefinitionReader(definitionFile, null);
    NamedAlleleMatcher namedAlleleMatcher = new NamedAlleleMatcher(new Env(), definitionReader, true, true, false);
    Result result = namedAlleleMatcher.call(new VcfFile(vcfFile), null);
    // ignore novel bases
    //printWarnings(result);
    assertEquals(3, result.getVcfWarnings().size());
    assertEquals(1, result.getGeneCalls().size());

    GeneCall geneCall = result.getGeneCalls().get(0);
    printMatches(geneCall);
    assertEquals(1, geneCall.getDiplotypes().size());
    DiplotypeMatch dm = geneCall.getDiplotypes().iterator().next();
    assertEquals("*2/[*17 + g.94781859G>A]", dm.getName());
    assertFalse(dm.getHaplotype1().getHaplotype().isCombination());
    assertFalse(dm.getHaplotype1().getHaplotype().isPartial());
    assertNotNull(dm.getHaplotype2());
    assertFalse(dm.getHaplotype2().getHaplotype().isCombination());
    assertTrue(dm.getHaplotype2().getHaplotype().isPartial());
  }

  @Test
  void testPartial3() throws Exception {
    Path definitionFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/NamedAlleleMatcher-partial3.json");
    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/NamedAlleleMatcher-partial3Phased.vcf");

    DefinitionReader definitionReader = new DefinitionReader(definitionFile, null);
    NamedAlleleMatcher namedAlleleMatcher = new NamedAlleleMatcher(new Env(), definitionReader, true, true, false);
    Result result = namedAlleleMatcher.call(new VcfFile(vcfFile), null);
    // ignore novel bases
    //printWarnings(result);
    assertEquals(6, result.getVcfWarnings().size());
    assertEquals(1, result.getGeneCalls().size());

    GeneCall geneCall = result.getGeneCalls().get(0);
    printMatches(geneCall);
    assertEquals(1, geneCall.getDiplotypes().size());
    DiplotypeMatch dm = geneCall.getDiplotypes().iterator().next();
    assertEquals("*10/[*9 + *14]", dm.getName());
    assertFalse(dm.getHaplotype1().getHaplotype().isCombination());
    assertFalse(dm.getHaplotype1().getHaplotype().isPartial());
    assertNotNull(dm.getHaplotype2());
    assertTrue(dm.getHaplotype2().getHaplotype().isCombination());
    assertFalse(dm.getHaplotype2().getHaplotype().isPartial());
  }

  @Test
  void testCombinationLongestScore() throws Exception {
    Path definitionFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/NamedAlleleMatcher-cyp2b6.json");
    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/NamedAlleleMatcher-combinationLongestScore.vcf");

    DefinitionReader definitionReader = new DefinitionReader(definitionFile, null);
    NamedAlleleMatcher namedAlleleMatcher = new NamedAlleleMatcher(new Env(), definitionReader, true, true, false);
    Result result = namedAlleleMatcher.call(new VcfFile(vcfFile), null);
    assertEquals(0, result.getVcfWarnings().size());
    assertEquals(1, result.getGeneCalls().size());

    GeneCall geneCall = result.getGeneCalls().get(0);
    List<String> matches = printMatches(geneCall);
    assertEquals(1, matches.size());
    assertThat(matches, contains("*13/[*6 + *14]"));

    DiplotypeMatch dm = geneCall.getDiplotypes().iterator().next();
    assertFalse(dm.getHaplotype1().getHaplotype().isCombination());
    assertFalse(dm.getHaplotype1().getHaplotype().isPartial());
    assertNotNull(dm.getHaplotype2());
    assertTrue(dm.getHaplotype2().getHaplotype().isCombination());
    assertFalse(dm.getHaplotype2().getHaplotype().isPartial());
  }


  /**
   * Make sure longest combination + partial scoring works.
   * {@code *1/*4 + *9 + g.41010006G>C} should beat {@code *1/*6 + g.41010006G>C}
   */
  @Test
  void testPartialLongestScore() throws Exception {
    Path definitionFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/NamedAlleleMatcher-cyp2b6.json");
    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/NamedAlleleMatcher-partialLongestScore.vcf");

    DefinitionReader definitionReader = new DefinitionReader(definitionFile, null);
    NamedAlleleMatcher namedAlleleMatcher = new NamedAlleleMatcher(new Env(), definitionReader, true, true, false);
    Result result = namedAlleleMatcher.call(new VcfFile(vcfFile), null);
    assertEquals(0, result.getVcfWarnings().size());
    assertEquals(1, result.getGeneCalls().size());

    GeneCall geneCall = result.getGeneCalls().get(0);
    List<String> matches = printMatches(geneCall);
    assertEquals(2, matches.size());
    assertThat(matches, contains("*1/[*6 + g.41010006G>C]", "*1/[*36 + g.41010006G>C]"));
  }

  @Test
  void testCombinationLongestScoreMissing(TestInfo testInfo) throws Exception {
    Path definitionFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/NamedAlleleMatcher-cyp2b6.json");
    Path vcfFile = new TestVcfBuilder(testInfo, "foo")
        .withDefinition(definitionFile)
        // c.1218G>A
        .variation("CYP2B6", "rs35773040", "G", "A")
        // c.1627A>G (*5)
        .variation("CYP2B6", "rs3745274", "T", "T")
        // c.85T>C (*9A)
        .variation("CYP2B6", "rs2279343", "G", "G")
        .missing("CYP2B6", "rs12721655")
        .generate();

    DefinitionReader definitionReader = new DefinitionReader(definitionFile, null);
    NamedAlleleMatcher namedAlleleMatcher = new NamedAlleleMatcher(new Env(), definitionReader, true, true, false);
    Result result = namedAlleleMatcher.call(new VcfFile(vcfFile), null);
    assertEquals(0, result.getVcfWarnings().size());
    assertEquals(1, result.getGeneCalls().size());

    GeneCall geneCall = result.getGeneCalls().get(0);
    printMatches(geneCall);
    assertEquals(4, geneCall.getDiplotypes().size());
    Optional<DiplotypeMatch> opt = geneCall.getDiplotypes().stream()
        .filter(dm -> dm.getName().equals("*13/[*6 + *14]"))
        .findFirst();
    assertTrue(opt.isPresent());
  }

  @Test
  void testPartialReferenceUnphased() throws Exception {
    Path definitionFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/NamedAlleleMatcher-cyp2b6.json");
    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/NamedAlleleMatcher-partialReferenceUnphased.vcf");

    DefinitionReader definitionReader = new DefinitionReader(definitionFile, null);
    NamedAlleleMatcher namedAlleleMatcher = new NamedAlleleMatcher(new Env(), definitionReader, true, true, false);
    Result result = namedAlleleMatcher.call(new VcfFile(vcfFile), null);
    // ignore novel bases
    //printWarnings(result);
    assertEquals(1, result.getVcfWarnings().size());
    assertEquals(1, result.getGeneCalls().size());

    GeneCall geneCall = result.getGeneCalls().get(0);
    printMatches(geneCall);
    assertEquals(1, geneCall.getDiplotypes().size());
    DiplotypeMatch dm = geneCall.getDiplotypes().iterator().next();
    assertEquals("*1/[*9 + g.41012316T>G]", dm.getName());
    assertFalse(dm.getHaplotype1().getHaplotype().isCombination());
    assertFalse(dm.getHaplotype1().getHaplotype().isPartial());
    assertNotNull(dm.getHaplotype2());
    assertFalse(dm.getHaplotype2().getHaplotype().isCombination());
    assertTrue(dm.getHaplotype2().getHaplotype().isPartial());
  }

  @Test
  void testPartialReferencePhased() throws Exception {
    Path definitionFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/NamedAlleleMatcher-cyp2b6.json");
    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/NamedAlleleMatcher-partialReferencePhased.vcf");

    DefinitionReader definitionReader = new DefinitionReader(definitionFile, null);
    NamedAlleleMatcher namedAlleleMatcher = new NamedAlleleMatcher(new Env(), definitionReader, true, true, false);
    Result result = namedAlleleMatcher.call(new VcfFile(vcfFile), null);
    // ignore novel bases
    //printWarnings(result);
    assertEquals(2, result.getVcfWarnings().size());
    assertEquals(1, result.getGeneCalls().size());

    GeneCall geneCall = result.getGeneCalls().get(0);
    printMatches(geneCall);
    assertEquals(1, geneCall.getDiplotypes().size());
    DiplotypeMatch dm = geneCall.getDiplotypes().iterator().next();
    assertEquals("*9/[g.41006968T>A + g.41012316T>G]", dm.getName());
    assertFalse(dm.getHaplotype1().getHaplotype().isCombination());
    assertFalse(dm.getHaplotype1().getHaplotype().isPartial());
    assertNotNull(dm.getHaplotype2());
    assertFalse(dm.getHaplotype2().getHaplotype().isCombination());
    assertTrue(dm.getHaplotype2().getHaplotype().isPartial());
  }


  @Test
  void testPartialReferenceDouble() throws Exception {
    Path definitionFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/NamedAlleleMatcher-cyp2b6.json");
    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/NamedAlleleMatcher-partialReferenceDouble.vcf");

    DefinitionReader definitionReader = new DefinitionReader(definitionFile, null);
    NamedAlleleMatcher namedAlleleMatcher = new NamedAlleleMatcher(new Env(), definitionReader, true, true, false);
    Result result = namedAlleleMatcher.call(new VcfFile(vcfFile), null);
    // ignore novel bases
    //printWarnings(result);
    assertEquals(1, result.getVcfWarnings().size());
    assertEquals(1, result.getGeneCalls().size());

    GeneCall geneCall = result.getGeneCalls().get(0);
    printMatches(geneCall);
    assertEquals(1, geneCall.getDiplotypes().size());
    DiplotypeMatch dm = geneCall.getDiplotypes().iterator().next();
    assertEquals("g.41006968T>A/g.41006968T>A", dm.getName());
    assertFalse(dm.getHaplotype1().getHaplotype().isCombination());
    assertTrue(dm.getHaplotype1().getHaplotype().isPartial());
    assertNotNull(dm.getHaplotype2());
    assertFalse(dm.getHaplotype2().getHaplotype().isCombination());
    assertTrue(dm.getHaplotype2().getHaplotype().isPartial());
  }


  @Test
  void testDpydEffectivelyPhased(TestInfo testInfo) throws Exception {
    Path definitionFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/NamedAlleleMatcher-dpyd.json");
    Path vcfFile = new TestVcfBuilder(testInfo, "Reference/c.62G>A")
        .withDefinition(definitionFile)
        // c.62G>A
        .variation("DPYD", "rs80081766", "C", "T")
        .generate();

    DefinitionReader definitionReader = new DefinitionReader(definitionFile, null);
    NamedAlleleMatcher namedAlleleMatcher = new NamedAlleleMatcher(new Env(), definitionReader, true, true, false);
    Result result = namedAlleleMatcher.call(new VcfFile(vcfFile), null);
    // ignore novel bases
    //printWarnings(result);
    assertEquals(0, result.getVcfWarnings().size());
    assertEquals(1, result.getGeneCalls().size());

    GeneCall geneCall = result.getGeneCalls().get(0);
    printMatches(geneCall);
    assertEquals(1, geneCall.getDiplotypes().size());
    DiplotypeMatch dm = geneCall.getDiplotypes().iterator().next();
    assertEquals("Reference/c.62G>A", dm.getName());
    System.out.println(geneCall.getHaplotypes());
    assertEquals(2, geneCall.getHaplotypes().size());
    List<String> names = geneCall.getHaplotypes().stream()
        .map(BaseMatch::getName)
        .toList();
    assertThat(names, contains("Reference", "c.62G>A"));
  }


  @Test
  void testDpydEffectivelyPhased2(TestInfo testInfo) throws Exception {
    // c.1218G>A/c.1218G>A c.1627A>G (*5)/c.1627A>G (*5) c.85T>C (*9A)
    Path definitionFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/NamedAlleleMatcher-dpyd.json");
    Path vcfFile = new TestVcfBuilder(testInfo, "[c.85T>C (*9A) + c.1218G>A + c.1627A>G (*5)]/[c.1218G>A + c.1627A>G (*5)]")
        .withDefinition(definitionFile)
        // c.1218G>A
        .variation("DPYD", "rs61622928", "T", "T")
        // c.1627A>G (*5)
        .variation("DPYD", "rs1801159", "C", "C")
        // c.85T>C (*9A)
        .variation("DPYD", "rs1801265", "A", "G")
        .generate();

    DefinitionReader definitionReader = new DefinitionReader(definitionFile, null);
    NamedAlleleMatcher namedAlleleMatcher = new NamedAlleleMatcher(new Env(), definitionReader, true, true, false);
    Result result = namedAlleleMatcher.call(new VcfFile(vcfFile), null);
    // ignore novel bases
    printWarnings(result);
    assertEquals(0, result.getVcfWarnings().size());
    assertEquals(1, result.getGeneCalls().size());

    GeneCall geneCall = result.getGeneCalls().get(0);
    printMatches(geneCall);
    assertEquals(1, geneCall.getDiplotypes().size());
    DiplotypeMatch dm = geneCall.getDiplotypes().iterator().next();
    assertEquals("[c.85T>C (*9A) + c.1218G>A + c.1627A>G (*5)]/[c.1218G>A + c.1627A>G (*5)]", dm.getName());
    System.out.println(geneCall.getHaplotypes());
    assertEquals(2, geneCall.getHaplotypes().size());
    Set<String> names = geneCall.getHaplotypes().stream()
        .flatMap(h -> {
          if (h instanceof CombinationMatch cm) {
            return cm.getComponentHaplotypes().stream()
                .map(NamedAllele::getName);
          }
          return Stream.of(h.getHaplotype().getName());
        })
        .collect(Collectors.toCollection(TreeSet::new));
    assertThat(new ArrayList<>(names), contains("c.1218G>A", "c.1627A>G (*5)", "c.85T>C (*9A)"));
  }


  @Test
  void testDpydPhased(TestInfo testInfo) throws Exception {
    Path definitionFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/NamedAlleleMatcher-dpyd.json");
    Path vcfFile = new TestVcfBuilder(testInfo, "Reference/c.62G>A + c.1129-5923C>G, c.1236G>A (HapB3) + c.3067C>A")
        .withDefinition(definitionFile)
        // c.62G>A
        .variation("DPYD", "rs80081766", "C", "T")
        // c.3067C>A
        .variation("DPYD", "rs114096998", "G", "T")
        .phased()
        .generate();

    DefinitionReader definitionReader = new DefinitionReader(definitionFile, null);
    NamedAlleleMatcher namedAlleleMatcher = new NamedAlleleMatcher(new Env(), definitionReader, true, true, false);
    Result result = namedAlleleMatcher.call(new VcfFile(vcfFile), null);
    // ignore novel bases
    //printWarnings(result);
    assertEquals(0, result.getVcfWarnings().size());
    assertEquals(1, result.getGeneCalls().size());

    GeneCall geneCall = result.getGeneCalls().get(0);
    printMatches(geneCall);
    assertEquals(1, geneCall.getDiplotypes().size());
    DiplotypeMatch dm = geneCall.getDiplotypes().iterator().next();
    assertEquals("Reference/[c.62G>A + c.3067C>A]", dm.getName());
    System.out.println(geneCall.getHaplotypes());
    assertEquals(2, geneCall.getHaplotypes().size());
    List<String> names = geneCall.getHaplotypes().stream()
        .map(BaseMatch::getName)
        .toList();
    assertThat(names, contains("Reference", "[c.62G>A + c.3067C>A]"));
  }

  @Test
  void testDpydUnPhased(TestInfo testInfo) throws Exception {
    Path definitionFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/NamedAlleleMatcher-dpyd.json");
    Path vcfFile = new TestVcfBuilder(testInfo, "Reference/c.62G>A + c.1129-5923C>G, c.1236G>A (HapB3) + c.3067C>A")
        .withDefinition(definitionFile)
        // c.62G>A
        .variation("DPYD", "rs80081766", "C", "T")
        // c.3067C>A
        .variation("DPYD", "rs114096998", "G", "T")
        .generate();

    DefinitionReader definitionReader = new DefinitionReader(definitionFile, null);
    NamedAlleleMatcher namedAlleleMatcher = new NamedAlleleMatcher(new Env(), definitionReader, true, true, false);
    Result result = namedAlleleMatcher.call(new VcfFile(vcfFile), null);
    // ignore novel bases
    //printWarnings(result);
    assertEquals(0, result.getVcfWarnings().size());
    assertEquals(1, result.getGeneCalls().size());

    GeneCall geneCall = result.getGeneCalls().get(0);
    printMatches(geneCall);
    assertEquals(0, geneCall.getDiplotypes().size());
    assertEquals(0, geneCall.getHaplotypes().size());
    assertEquals(2, geneCall.getHaplotypeMatches().size());
    List<String> names = geneCall.getHaplotypeMatches().stream()
        .map(BaseMatch::getName)
        .toList();
    assertThat(names, contains("c.62G>A", "c.3067C>A"));
  }

  @Test
  void testDpydPhasedDouble(TestInfo testInfo) throws Exception {
    Path definitionFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/NamedAlleleMatcher-dpyd.json");
    Path vcfFile = new TestVcfBuilder(testInfo, "c.62G>A + c.3067C>A/c.62G>A + c.1129-5923C>G, c.1236G>A (HapB3) + c.3067C>A")
        .withDefinition(definitionFile)
        // c.1129-5923C>G, c.1236G>A (HapB3)
        .variation("DPYD", "rs75017182", "G", "C")
        .variation("DPYD", "rs56038477", "C", "T")
        // c.62G>A
        .variation("DPYD", "rs80081766", "T", "T")
        // c.3067C>A
        .variation("DPYD", "rs114096998", "T", "T")
        .phased()
        .generate();

    DefinitionReader definitionReader = new DefinitionReader(definitionFile, null);
    NamedAlleleMatcher namedAlleleMatcher = new NamedAlleleMatcher(new Env(), definitionReader, true, true, false);
    Result result = namedAlleleMatcher.call(new VcfFile(vcfFile), null);
    // ignore novel bases
    //printWarnings(result);
    assertEquals(0, result.getVcfWarnings().size());
    assertEquals(1, result.getGeneCalls().size());

    GeneCall geneCall = result.getGeneCalls().get(0);
    System.out.println("diplotypes:");
    printMatches(geneCall);
    System.out.println("\nhaplotypes:");
    System.out.println(geneCall.getHaplotypes());
    assertEquals(1, geneCall.getDiplotypes().size());
    DiplotypeMatch dm = geneCall.getDiplotypes().iterator().next();
    assertEquals("[c.62G>A + c.1129-5923C>G, c.1236G>A (HapB3) + c.3067C>A]/[c.62G>A + c.3067C>A]", dm.getName());
    assertEquals(2, geneCall.getHaplotypes().size());
    List<String> names = geneCall.getHaplotypes().stream()
        .map(BaseMatch::getName)
        .toList();
    assertThat(names, contains("[c.62G>A + c.1129-5923C>G, c.1236G>A (HapB3) + c.3067C>A]", "[c.62G>A + c.3067C>A]"));
  }


  @Test
  void testDpydHomozygous(TestInfo testInfo) throws Exception {
    Path definitionFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/NamedAlleleMatcher-dpyd.json");
    Path vcfFile = new TestVcfBuilder(testInfo, "c.62G>A + c.3067C>A/c.62G>A + c.3067C>A")
        .withDefinition(definitionFile)
        // c.62G>A
        .variation("DPYD", "rs80081766", "T", "T")
        // c.3067C>A
        .variation("DPYD", "rs114096998", "T", "T")
        .phased()
        .generate();

    DefinitionReader definitionReader = new DefinitionReader(definitionFile, null);
    NamedAlleleMatcher namedAlleleMatcher = new NamedAlleleMatcher(new Env(), definitionReader, true, true, false);
    Result result = namedAlleleMatcher.call(new VcfFile(vcfFile), null);
    // ignore novel bases
    //printWarnings(result);
    assertEquals(0, result.getVcfWarnings().size());
    assertEquals(1, result.getGeneCalls().size());

    GeneCall geneCall = result.getGeneCalls().get(0);
    printMatches(geneCall);
    System.out.println(geneCall.getHaplotypes());
    assertEquals(1, geneCall.getDiplotypes().size());
    DiplotypeMatch dm = geneCall.getDiplotypes().iterator().next();
    assertEquals("[c.62G>A + c.3067C>A]/[c.62G>A + c.3067C>A]", dm.getName());
    List<String> names = geneCall.getHaplotypes().stream()
        .map(BaseMatch::getName)
        .toList();
    assertThat(names, contains("[c.62G>A + c.3067C>A]"));
  }


  @Test
  void testDpydEffectivelyPhasedCombination(TestInfo testInfo) throws Exception {
    // this test is based on NA18973
    Path definitionFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/NamedAlleleMatcher-dpyd.json");
    Path vcfFile = new TestVcfBuilder(testInfo, "c.1627A>G (*5)/c.1627A>G (*5) + c.85T>C (*9A)")
        .withDefinition(definitionFile)
        .variation("DPYD", "rs1801159", "C", "C")
        .variation("DPYD", "rs1801265", "A", "G")
        .missing("DPYD", "rs148799944", "rs140114515", "rs1801268", "rs72547601", "rs72547602",
            "rs141044036", "rs147545709", "rs55674432", "rs146529561", "rs137999090", "rs138545885", "rs55971861",
            "rs72549303", "rs147601618", "rs145773863", "rs138616379", "rs148994843", "rs138391898", "rs111858276",
            "rs72549304", "rs142512579", "rs764666241", "rs140602333", "rs78060119", "rs143154602", "rs72549306",
            "rs145112791", "rs150437414", "rs1801266", "rs72549307", "rs72549308", "rs139834141", "rs141462178",
            "rs150385342", "rs72549309", "rs80081766", "rs150036960")
        .generate();

    DefinitionReader definitionReader = new DefinitionReader(definitionFile, null);
    NamedAlleleMatcher namedAlleleMatcher = new NamedAlleleMatcher(new Env(), definitionReader, true, true, false);
    Result result = namedAlleleMatcher.call(new VcfFile(vcfFile), null);
    // ignore novel bases
    //printWarnings(result);
    assertEquals(0, result.getVcfWarnings().size());
    assertEquals(1, result.getGeneCalls().size());

    GeneCall geneCall = result.getGeneCalls().get(0);
    printMatches(geneCall);
    System.out.println(geneCall.getHaplotypes());
    assertEquals(1, geneCall.getDiplotypes().size());
    DiplotypeMatch dm = geneCall.getDiplotypes().iterator().next();
    assertEquals("c.1627A>G (*5)/[c.85T>C (*9A) + c.1627A>G (*5)]", dm.getName());
    assertEquals(2, geneCall.getHaplotypes().size());
    List<String> names = geneCall.getHaplotypes().stream()
        .map(BaseMatch::getName)
        .toList();
    assertThat(names, contains("c.1627A>G (*5)", "[c.85T>C (*9A) + c.1627A>G (*5)]"));
  }


  @Test
  void testDpydUnphasedHomozygousNoFunction(TestInfo testInfo) throws Exception {
    Path definitionFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/NamedAlleleMatcher-dpyd.json");
    Path vcfFile = new TestVcfBuilder(testInfo, "c.1627A>G (*5)/c.1627A>G (*5) + c.85T>C (*9A)")
        .withDefinition(definitionFile)
        .variation("DPYD", "rs72547601", "C", "C") // c.2933A>G
        .variation("DPYD", "rs67376798", "A", "T") // c.2846A>T
        .variation("DPYD", "rs60139309", "T", "C") // c.2582A>G
        .generate();

    DefinitionReader definitionReader = new DefinitionReader(definitionFile, null);
    NamedAlleleMatcher namedAlleleMatcher = new NamedAlleleMatcher(new Env(), definitionReader, true, true, false);
    Result result = namedAlleleMatcher.call(new VcfFile(vcfFile), null);
    // ignore novel bases
    //printWarnings(result);
    assertEquals(0, result.getVcfWarnings().size());
    assertEquals(1, result.getGeneCalls().size());

    GeneCall geneCall = result.getGeneCalls().get(0);
    printMatches(geneCall);
    System.out.println(geneCall.getHaplotypes());
    assertEquals(0, geneCall.getDiplotypes().size());
    assertEquals(0, geneCall.getHaplotypes().size());
    assertEquals(4, geneCall.getHaplotypeMatches().size());
    assertThat(geneCall.getHaplotypeMatches().stream().map(BaseMatch::toString).toList(),
        contains("c.2582A>G", "c.2846A>T", "c.2933A>G", "c.2933A>G"));
  }


  /**
   * This tests that sorting {@link DiplotypeMatch} works correctly.
   * Problem originally found by Andrew while testing DPYD.
   */
  @Test
  void testDiplotypeMatcher() throws Exception {
    Path definitionFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/NamedAlleleMatcher-dpyd.json");
    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/NamedAlleleMatcher-diplotypeMatcher.vcf");

    DefinitionReader definitionReader = new DefinitionReader(definitionFile, null);
    // this problem doesn't happen consistently, which is why we are doing this in a loop
    for (int x = 0; x < 10; x += 1) {
      NamedAlleleMatcher namedAlleleMatcher = new NamedAlleleMatcher(new Env(), definitionReader, true, true, false);
      Result result = namedAlleleMatcher.call(new VcfFile(vcfFile), null);
      assertEquals(1, result.getGeneCalls().size());
    }
  }


  @Test
  void testUnknownAltMultisample() throws Exception {
    Path definitionFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/NamedAlleleMatcher-ryr1.json");
    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/NamedAlleleMatcher-unknownAltMultisample.vcf");

    DefinitionReader definitionReader = new DefinitionReader(definitionFile, null);
    NamedAlleleMatcher namedAlleleMatcher = new NamedAlleleMatcher(new Env(), definitionReader, true, true, false);
    Result result = namedAlleleMatcher.call(new VcfFile(vcfFile), "Sample_2");
    // ignore novel bases
    printWarnings(result);
    assertEquals(0, result.getVcfWarnings().size());
  }


  @SuppressWarnings("unused")
  private static void printWarnings(Result result) {
    for (String key : result.getVcfWarnings().keySet()) {
      System.out.println(key);
      System.out.println("\t" + result.getVcfWarnings().get(key));
    }
  }

  private static List<String> printMatches(GeneCall geneCall) {
    List<String> names = new ArrayList<>();
    for (DiplotypeMatch d : geneCall.getDiplotypes()) {
      System.out.println(d.getName() + " (" + d.getScore() + ": " + d.getHaplotype1().getHaplotype().getScore() +
          " / " + Objects.requireNonNull(d.getHaplotype2()).getHaplotype().getScore() + ")");
      names.add(d.getName());
    }
    return names;
  }
}

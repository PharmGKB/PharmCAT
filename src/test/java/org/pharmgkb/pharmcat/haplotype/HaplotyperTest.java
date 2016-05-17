package org.pharmgkb.pharmcat.haplotype;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import javax.annotation.Nonnull;
import com.google.common.base.Joiner;
import com.google.common.base.Preconditions;
import com.google.common.collect.Sets;
import org.junit.Test;
import org.pharmgkb.common.util.PathUtils;
import org.pharmgkb.pharmcat.TestUtil;
import org.pharmgkb.pharmcat.definition.model.NamedAllele;
import org.pharmgkb.pharmcat.haplotype.model.DiplotypeMatch;
import org.pharmgkb.pharmcat.haplotype.model.HaplotyperResult;

import static org.hamcrest.core.IsEqual.equalTo;
import static org.junit.Assert.*;


/**
 * JUnit test for {@link Haplotyper}.
 *
 * @author Mark Woon
 */
public class HaplotyperTest {


  static HaplotyperResult testCallHaplotype(@Nonnull Path tsvFile, @Nonnull Path vcfFile) throws Exception {
    return testCallHaplotype(tsvFile, vcfFile, true, true, false);
  }

  /**
   * Helper method for running Haplotyper.
   * This is used by the more specific gene tests.
   */
  static HaplotyperResult testCallHaplotype(@Nonnull Path definitionFile, @Nonnull Path vcfFile,
      boolean assumeReference, boolean topCandidateOnly, boolean showUnmatched) throws Exception {

    DefinitionReader definitionReader = new DefinitionReader();
    definitionReader.read(definitionFile);

    Haplotyper haplotyper = new Haplotyper(definitionReader, assumeReference, topCandidateOnly);
    HaplotyperResult result = haplotyper.call(vcfFile);

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
   * @param expectedPairs the set of expected diplotypes in "*1/*2" format
   * @param result the {@link Haplotyper} results
   */
  static void assertDiplotypePairs(@Nonnull List<String> expectedPairs,
      @Nonnull HaplotyperResult result) {

    Preconditions.checkNotNull(expectedPairs);
    Preconditions.checkNotNull(result);

    List<String> pairs = new ArrayList<>();
    for (DiplotypeMatch dm : result.getGeneCalls().get(0).getDiplotypes()) {
      pairs.add(dm.getName());
    }

    if (expectedPairs.size() != pairs.size() || !expectedPairs.equals(pairs)) {
      System.out.println("Expected: [" + Joiner.on(", ").join(expectedPairs));
      System.out.println("Got:      " + pairs);
      fail("Did not get expected matches");
    }
  }


  @Test
  public void testCall() throws Exception {

    Path vcfFile  = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/haplotyper.vcf");
    Path jsonFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/haplotyper.json");

    DefinitionReader definitionReader = new DefinitionReader();
    definitionReader.read(jsonFile);

    Haplotyper haplotyper = new Haplotyper(definitionReader);
    HaplotyperResult result = haplotyper.call(vcfFile);
    Set<DiplotypeMatch> pairs = result.getGeneCalls().get(0).getDiplotypes();
    assertNotNull(pairs);
    assertEquals(1, pairs.size());
    assertEquals("*1/*2", pairs.iterator().next().getName());
  }


  /**
   * This breaks down the main code path that {@link #testCall()} runs to simplify testing smaller chunks at a time.
   */
  @Test
  public void testCallDiplotypePath() throws Exception {

    Path vcfFile  = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/haplotyper.vcf");
    Path jsonFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/haplotyper.json");
    Set<String> permutations = Sets.newHashSet(
        "1:C;2:del;3:C;",
        "1:T;2:insA;3:delC;"
    );
    String gene = "CYP3A5";

    DefinitionReader definitionReader = new DefinitionReader();
    definitionReader.read(jsonFile);

    Haplotyper haplotyper = new Haplotyper(definitionReader);
    VcfReader vcfReader = haplotyper.buildVcfReader(vcfFile);

    // grab SampleAlleles for all positions related to current gene
    MatchData data = new MatchData(vcfReader.getAlleleMap(), "chr1", definitionReader.getPositions(gene));
    assertEquals(3, data.getNumSampleAlleles());
    assertEquals(0, data.getMissingPositions().size());
    // handle missing positions of interest in sample
    data.marshallHaplotypes(definitionReader.getHaplotypes(gene));
    assertEquals(3, data.getPositions().length);
    assertEquals(2, data.getHaplotypes().size());

    // get all permutations of sample at positions of interest
    data.generateSamplePermutations();
    assertThat(data.getPermutations(), equalTo(permutations));

    for (NamedAllele hap : data.getHaplotypes()) {
      System.out.println(hap.getName() + ": " + hap.getPermutations().pattern());
    }

    List<DiplotypeMatch> pairs = new DiplotypeMatcher(data).compute();
    assertNotNull(pairs);
    assertEquals(1, pairs.size());
    assertEquals("*1/*2", pairs.get(0).getName());
  }
}

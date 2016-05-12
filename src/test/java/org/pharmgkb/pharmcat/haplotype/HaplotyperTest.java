package org.pharmgkb.pharmcat.haplotype;

import java.nio.file.Path;
import java.util.List;
import java.util.Set;
import java.util.SortedMap;
import javax.annotation.Nonnull;
import org.junit.Test;
import org.pharmgkb.common.util.PathUtils;
import org.pharmgkb.pharmcat.TestUtil;
import org.pharmgkb.pharmcat.definition.model.NamedAllele;
import org.pharmgkb.pharmcat.haplotype.model.DiplotypeMatch;
import org.pharmgkb.pharmcat.haplotype.model.HaplotyperResult;

import static org.junit.Assert.*;


/**
 * JUnit test for {@link Haplotyper}.
 *
 * @author Mark Woon
 */
public class HaplotyperTest {


  static List<DiplotypeMatch> testCallHaplotype(@Nonnull Path tsvFile, @Nonnull Path vcfFile) throws Exception {
    return testCallHaplotype(tsvFile, vcfFile, false);
  }

  /**
   * Helper method for running Haplotyper.
   * This is used by the more specific gene tests.
   */
  static List<DiplotypeMatch> testCallHaplotype(@Nonnull Path tsvFile, @Nonnull Path vcfFile, boolean showUnmatched)
      throws Exception {

    DefinitionReader definitionReader = new DefinitionReader();
    definitionReader.read(tsvFile);
    String gene = definitionReader.getGenes().iterator().next();

    Haplotyper haplotyper = new Haplotyper(definitionReader);
    SortedMap<String, SampleAllele> alleleMap = haplotyper.getVcfReader().read(vcfFile);

    MatchData dataset = haplotyper.initializeCallData(alleleMap, gene);
    List<DiplotypeMatch> matches = haplotyper.callDiplotypes(dataset);
    StringBuilder rezBuilder = new StringBuilder();
    for (DiplotypeMatch dm : matches) {
      if (rezBuilder.length() > 0) {
        rezBuilder.append(", ");
      }
      rezBuilder.append(dm.getName())
          .append(" (")
          .append(dm.getScore())
          .append(")");
    }
    System.out.println(rezBuilder);

    HaplotyperResult result = new ResultBuilder(definitionReader)
        .forFile(vcfFile)
        .gene(gene, dataset, matches)
        .build();
    // print
    new ResultSerializer()
        .alwaysShowUnmatchedHaplotypes(showUnmatched)
        .toHtml(result, vcfFile.getParent().resolve(PathUtils.getBaseFilename(vcfFile) + ".html"))
        .toJson(result, vcfFile.getParent().resolve(PathUtils.getBaseFilename(vcfFile) + ".json"));

    return matches;
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

    String gene = "CYP3A5";

    DefinitionReader definitionReader = new DefinitionReader();
    definitionReader.read(jsonFile);

    Haplotyper haplotyper = new Haplotyper(definitionReader);
    SortedMap<String, SampleAllele> alleles = haplotyper.getVcfReader().read(vcfFile);

    MatchData data = new MatchData();
    // grab SampleAlleles for all positions related to current gene
    data.marshallSampleData(alleles, "chr1", definitionReader.getPositions(gene));
    assertEquals(3, data.geneSampleMap.size());
    assertEquals(0, data.missingPositions.size());
    // handle missing positions of interest in sample
    data.marshallHaplotypes(definitionReader.getPositions(gene), definitionReader.getHaplotypes(gene));
    assertEquals(3, data.positions.length);
    assertEquals(2, data.haplotypes.size());

    // get all permutations of sample at positions of interest
    data.generateSamplePermutations();
    assertEquals(2, data.permutations.size());
    assertTrue(data.permutations.contains("1:C;2:del;3:C;"));
    assertTrue(data.permutations.contains("1:T;2:insA;3:delC;"));

    for (NamedAllele hap : data.haplotypes) {
      System.out.println(hap.getName() + ": " + hap.getPermutations().pattern());
    }

    List<DiplotypeMatch> pairs = new DiplotypeMatcher(data).compute();
    assertNotNull(pairs);
    assertEquals(1, pairs.size());
    assertEquals("*1/*2", pairs.get(0).getName());
  }
}

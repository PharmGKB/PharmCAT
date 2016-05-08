package org.pharmgkb.pharmcat.haplotype;

import java.nio.file.Path;
import java.util.List;
import java.util.Set;
import java.util.SortedMap;
import java.util.stream.Collectors;
import javax.annotation.Nonnull;
import org.junit.Test;
import org.pharmgkb.pharmcat.TestUtil;
import org.pharmgkb.pharmcat.definition.model.NamedAllele;
import org.pharmgkb.pharmcat.haplotype.model.DiplotypeMatch;

import static org.junit.Assert.*;


/**
 * JUnit test for {@link Haplotyper}.
 *
 * @author Mark Woon
 */
public class HaplotyperTest {

  /**
   * Helper method for running Haplotyper.
   * This is used by the more specific gene tests.
   */
  static List<DiplotypeMatch> testCallHaplotype(@Nonnull Path tsvFile, @Nonnull Path vcfFile) throws Exception {

    DefinitionReader definitionReader = new DefinitionReader();
    definitionReader.read(tsvFile);
    String gene = definitionReader.getGenes().iterator().next();

    Haplotyper haplotyper = new Haplotyper(definitionReader);
    SortedMap<String, SampleAllele> alleleMap = haplotyper.getVcfReader().read(vcfFile);

    List<DiplotypeMatch> matches = haplotyper.callDiplotypes(alleleMap, gene);
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

    // print
    new Report(definitionReader)
        .forFile(vcfFile)
        .gene(gene, matches, alleleMap.values())
        .printHtml();

    return matches;
  }


  @Test
  public void testCall() throws Exception {

    Path vcfFile  = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/haplotyper.vcf");
    Path jsonFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/haplotyper.json");

    DefinitionReader definitionReader = new DefinitionReader();
    definitionReader.read(jsonFile);

    Haplotyper haplotyper = new Haplotyper(definitionReader);
    Report report = haplotyper.call(vcfFile);
    Set<DiplotypeMatch> pairs = report.getResults().getGeneCalls().get(0).getDiplotypes();
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

    Haplotyper.Dataset data = new Haplotyper.Dataset();
    // grab SampleAlleles for all positions related to current gene
    data.marshallSampleData(alleles, "chr1", definitionReader.getPositions(gene));
    assertEquals(3, data.geneSampleMap.size());
    assertEquals(0, data.missingPositions.size());
    // handle missing positions of interest in sample
    data.marshallHaplotypes(definitionReader.getPositions(gene), definitionReader.getHaplotypes(gene));
    assertEquals(3, data.positions.length);
    assertEquals(2, data.haplotypes.size());

    // get all permutations of sample at positions of interest
    Set<String> permutations = CombinationUtil.generatePermutations(
        data.geneSampleMap.values().stream()
            .sorted()
            .collect(Collectors.toList())
    );
    assertEquals(2, permutations.size());
    assertTrue(permutations.contains("1:C;2:del;3:C;"));
    assertTrue(permutations.contains("1:T;2:insA;3:delC;"));

    for (NamedAllele hap : data.haplotypes) {
      System.out.println(hap.getName() + ": " + hap.getPermutations().pattern());
    }

    List<DiplotypeMatch> pairs = new DiplotypeMatcher(data.geneSampleMap, permutations, data.haplotypes).compute();
    assertNotNull(pairs);
    assertEquals(1, pairs.size());
    assertEquals("*1/*2", pairs.get(0).getName());
  }
}

package org.pharmgkb.pharmcat.haplotype;

import java.nio.file.Path;
import java.util.List;
import java.util.SortedMap;
import javax.annotation.Nonnull;
import org.junit.Test;
import org.pharmgkb.pharmcat.TestUtil;
import org.pharmgkb.pharmcat.haplotype.model.DiplotypeMatch;


/**
 * JUnit test for {@link Haplotyper}.
 *
 * @author Mark Woon
 */
public class HaplotyperTest {


  static List<DiplotypeMatch> testCallHaplotype(@Nonnull Path tsvFile, @Nonnull Path vcfFile) throws Exception {

    DefinitionReader definitionReader = new DefinitionReader();
    definitionReader.read(tsvFile);
    String gene = definitionReader.getHaplotypes().keySet().iterator().next();

    VcfReader vcfReader = new VcfReader(Haplotyper.calculateLocationsOfInterest(definitionReader));
    SortedMap<String, SampleAllele> alleleMap = vcfReader.read(vcfFile);

    Haplotyper haplotyper = new Haplotyper(definitionReader);
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
  public void testCyp2c19() throws Exception {
    // initial test - file reading etc
    Path tsvFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/CYP2C19.tsv");
    Path vcfFile = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/NA12878.2c19_filtered.vcf");

    DefinitionReader definitionReader = new DefinitionReader();
    definitionReader.read(tsvFile);

    Haplotyper haplotyper = new Haplotyper(definitionReader);
    Report report = haplotyper.call(vcfFile);
  }
}

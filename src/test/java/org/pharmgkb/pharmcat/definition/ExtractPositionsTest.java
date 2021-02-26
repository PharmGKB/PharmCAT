package org.pharmgkb.pharmcat.definition;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.List;
import com.google.common.collect.ImmutableList;
import org.junit.jupiter.api.Test;
import org.pharmgkb.common.util.PathUtils;
import org.pharmgkb.pharmcat.haplotype.DefinitionReader;

import static org.junit.jupiter.api.Assertions.assertEquals;

/**
 * JUnit test for {@link ExtractPositions}.
 *
 * @author lester
 */
class ExtractPositionsTest {

  /**
   * Test to see if making a VCF of just VKORC1 positions will have the same count of positions as the input
   * definition file.
   * @throws IOException can occur when writing temporary files to the filesystem
   */
  @Test
  void testExtract() throws IOException {
    List<String> genesToTest = ImmutableList.of("VKORC1", "CYP2C9");

    DefinitionReader definitionReader = new DefinitionReader();
    for (String geneToTest : genesToTest) {
      definitionReader.read(PathUtils.getPathToResource("org/pharmgkb/pharmcat/definition/alleles/" + geneToTest + "_translation.json").toAbsolutePath());
    }
    Path outputVcf = Files.createTempFile("positions", ".vcf");
    try {
      ExtractPositions extractPositions = new ExtractPositions(outputVcf);
      String[] vcfLines = extractPositions.getPositions().toString().split("\\n");

      for (String geneToTest : genesToTest) {
        int numPositionsInDefinitionFile = definitionReader.getPositions(geneToTest).length;
        long numPositionsInVcfFile = Arrays.stream(vcfLines).filter(l -> !l.startsWith("#") && l.contains(geneToTest)).count();
        assertEquals(numPositionsInDefinitionFile, numPositionsInVcfFile, geneToTest + " has mismatched number of positions between definition and VCF");
      }
    } finally {
      Files.deleteIfExists(outputVcf);
    }
  }

  @Test
  void testGetDAS() throws IOException {
    ExtractPositions extractPositions = new ExtractPositions(Files.createTempFile("testGetDAS", "vcf"));
    String ref = extractPositions.getDAS("chr6", "18149127", "hg38");
    assertEquals("T", ref);
  }
}

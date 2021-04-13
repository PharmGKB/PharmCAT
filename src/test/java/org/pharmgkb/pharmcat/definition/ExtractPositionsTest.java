package org.pharmgkb.pharmcat.definition;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Map;
import org.junit.jupiter.api.Test;
import org.pharmgkb.pharmcat.haplotype.DefinitionReader;
import org.pharmgkb.pharmcat.util.DataManager;

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
    DefinitionReader definitionReader = new DefinitionReader();
    definitionReader.read(DataManager.DEFAULT_DEFINITION_DIR);
    Path outputVcf = Files.createTempFile("positions", ".vcf");
    try {
      ExtractPositions extractPositions = new ExtractPositions(outputVcf);
      Map<Integer, Map<Integer, String[]>> chrMap = extractPositions.getPositions();

      for (String gene : definitionReader.getGenes()) {
        int numPositionsInDefinitionFile = definitionReader.getPositions(gene).length;
        long numPositionsInVcfFile = chrMap.values().stream()
            .flatMap(m -> m.values().stream())
            .map(s -> s[7])
            .filter(s -> s.contains(gene))
            .count();
        assertEquals(numPositionsInDefinitionFile, numPositionsInVcfFile, gene + " has " +
            numPositionsInVcfFile + " defined positions but only " + numPositionsInVcfFile + " can be found in VCF");
      }
    } finally {
      Files.deleteIfExists(outputVcf);
    }
  }

  @Test
  void testPreviousBase() throws IOException {
    ExtractPositions extractPositions = new ExtractPositions(Files.createTempFile("testGetDAS", "vcf"));
    String ref = extractPositions.getPreviousBase("hg38", "chr6", 18149127);
    assertEquals("T", ref);
  }
}

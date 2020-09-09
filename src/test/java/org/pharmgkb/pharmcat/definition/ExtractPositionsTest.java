package org.pharmgkb.pharmcat.definition;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
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
    DefinitionReader definitionReader = new DefinitionReader();
    File file = PathUtils.getPathToResource("org/pharmgkb/pharmcat/definition/alleles/VKORC1_translation.json").toFile();
    Path path = Paths.get(file.getAbsolutePath());
    definitionReader.read(path);
    Path outputVcf = Files.createTempFile("positions", ".vcf");
    try {
      ExtractPositions extractPositions = new ExtractPositions(outputVcf);
      StringBuilder vcfText = extractPositions.getPositions(definitionReader);

      int numPositionsInDefinitinoFile = definitionReader.getPositions("VKORC1").length;
      long numPositionsInVcfFile = Arrays.stream(vcfText.toString().split("\\n")).filter(l -> !l.startsWith("#")).count();
      assertEquals(numPositionsInDefinitinoFile, numPositionsInVcfFile);
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

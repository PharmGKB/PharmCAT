package org.pharmgkb.pharmcat.definition;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import org.apache.commons.lang3.StringUtils;
import org.junit.Test;
import org.pharmgkb.common.util.PathUtils;
import org.pharmgkb.pharmcat.haplotype.DefinitionReader;

import static org.junit.Assert.assertEquals;


/**
 *
 * JUnit test for {@link ExtractPositions}.
 *
 * Created by lester on 6/21/16.
 */
public class ExtractPositionsTest {


  // Quick test on a single definition file, comparing number of positions with final number of variants in vcf string
  @Test
  public void testExtract() throws IOException {
    DefinitionReader definitionReader = new DefinitionReader();
    File file = PathUtils.getPathToResource("org/pharmgkb/pharmcat/definition/alleles/VKORC1_translation.json").toFile();
    Path path = Paths.get(file.getAbsolutePath());
    definitionReader.read(path);
    Path outputVcf = Paths.get("org/pharmgkb/pharmcat/haplotype/positions.vcf");
    ExtractPositions extractPositions = new ExtractPositions(path, outputVcf, "hg38");
    StringBuilder vcfText = extractPositions.getPositions(definitionReader);
    int definitionPositions = definitionReader.getPositions("VKORC1").length;
    int vcfVariants = StringUtils.countMatches(vcfText.toString(), "chr");
    assertEquals(definitionPositions, vcfVariants);
  }


}

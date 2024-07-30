package org.pharmgkb.pharmcat.haplotype;

import java.nio.file.Path;
import org.junit.jupiter.api.Test;
import org.pharmgkb.common.util.PathUtils;
import org.pharmgkb.pharmcat.definition.DefinitionReader;
import org.pharmgkb.pharmcat.definition.model.VariantLocus;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;


class DefinitionReaderTest {


  @Test
  void testCyp3a4() throws Exception {
    Path definitionFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/definition/alleles/CYP3A4_translation.json");
    DefinitionReader dr = new DefinitionReader(definitionFile, null);

    assertEquals(2, dr.getPositions("CYP3A4").length);
    assertEquals(3, dr.getHaplotypes("CYP3A4").size());
  }


  @Test
  void testReadAllDefinitions() throws Exception {

    DefinitionReader reader = DefinitionReader.defaultReader();

    for (String gene : reader.getGenes()) {
      assertTrue(reader.getDefinitionFile(gene).getChromosome().startsWith("chr"));
      for (VariantLocus variant : reader.getPositions(gene)) {
        assertTrue(variant.getPosition() > 0, "Zero position [" + variant.getPosition() + "] for variant " + variant.getChromosomeHgvsName() + " in gene " + gene);
      }
    }
  }
}

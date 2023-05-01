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
  void testVkorc1() throws Exception {
    Path definitionFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/definition/alleles/VKORC1_translation.json");
    DefinitionReader dr = new DefinitionReader(definitionFile, null);

    assertEquals(1, dr.getPositions("VKORC1").length);
    assertEquals(2, dr.getHaplotypes("VKORC1").size());
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

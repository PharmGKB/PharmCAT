package org.pharmgkb.pharmcat.haplotype;

import java.io.File;
import java.nio.file.Path;
import java.nio.file.Paths;
import org.junit.jupiter.api.Test;
import org.pharmgkb.common.util.PathUtils;
import org.pharmgkb.pharmcat.definition.DefinitionReader;
import org.pharmgkb.pharmcat.definition.model.VariantLocus;
import org.pharmgkb.pharmcat.util.DataManager;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;


class DefinitionReaderTest {


  @Test
  void testVkorc1() throws Exception {
    DefinitionReader dr = new DefinitionReader();
    File file = PathUtils.getPathToResource("org/pharmgkb/pharmcat/definition/alleles/VKORC1_translation.json").toFile();
    Path path = Paths.get(file.getAbsolutePath());
    dr.read(path);

    assertEquals(1, dr.getPositions("VKORC1").length);
    assertEquals(2, dr.getHaplotypes("VKORC1").size());
  }


  @Test
  void testReadAllDefinitions() throws Exception {

    DefinitionReader reader = new DefinitionReader();
    reader.read(DataManager.DEFAULT_DEFINITION_DIR);

    for (String gene : reader.getGenes()) {
      assertTrue(reader.getDefinitionFile(gene).getChromosome().startsWith("chr"));
      for (VariantLocus variant : reader.getPositions(gene)) {
        assertTrue(variant.getPosition() > 0, "Zero position [" + variant.getPosition() + "] for variant " + variant.getChromosomeHgvsName() + " in gene " + gene);
      }
    }
  }
}

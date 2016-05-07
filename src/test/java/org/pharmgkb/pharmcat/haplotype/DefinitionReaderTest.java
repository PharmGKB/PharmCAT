package org.pharmgkb.pharmcat.haplotype;

import java.io.File;
import java.nio.file.Path;
import java.nio.file.Paths;
import org.junit.Test;
import org.pharmgkb.pharmcat.TestUtil;
import org.pharmgkb.pharmcat.definition.model.VariantLocus;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;


public class DefinitionReaderTest {


  @Test
  public void testVkorc1() throws Exception {

    System.out.println("DefinitionReaderTest");

    DefinitionReader dr = new DefinitionReader();
    File file = new File(DefinitionReader.class.getResource("VKORC1.json").getFile());
    Path path = Paths.get(file.getAbsolutePath());
    dr.read(path);

    assertEquals(1, dr.getPositions("VKORC1").length);
    assertEquals(2, dr.getHaplotypes("VKORC1").size());
  }


  @Test
  public void testReadAllDefinitions() throws Exception {

    Path file = TestUtil.getFile("org/pharmgkb/pharmcat/haplotype/CYP2C19.json");
    DefinitionReader reader = new DefinitionReader();
    reader.read(file.getParent());

    for (String gene : reader.getGenes()) {
      assertTrue(reader.getDefinitionFile(gene).getChromosome().startsWith("chr"));
      for (VariantLocus variant : reader.getPositions(gene)) {
        assertTrue(variant.getPosition() > 0);
      }
    }
  }
}

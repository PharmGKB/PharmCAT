package org.pharmgkb.pharmcat.util;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import org.junit.Test;
import org.pharmgkb.common.util.PathUtils;
import org.pharmgkb.pharmcat.definition.model.VariantLocus;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;


/**
 * Test for Ensembl utility class
 *
 * @author Ryan Whaley
 */
public class EnsemblUtilsTest {
  private static final Path sf_variantTestDataPath = PathUtils.getPathToResource("org/pharmgkb/pharmcat/util/EnsemblUtilsVariant.json");

  /**
   * Testing a known variant record
   */
  @Test
  public void testParse() throws IOException {
    VariantLocus locus = EnsemblUtils.parse(Files.newInputStream(sf_variantTestDataPath));

    assertNotNull(locus);
    assertEquals("rs12777823", locus.getRsid());
    assertEquals(94645745, locus.getPosition());
    assertEquals("chr10", locus.getChromosome());
    assertEquals("chr10:94645745", locus.getVcfChrPosition());
  }
}

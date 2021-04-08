package org.pharmgkb.pharmcat.haplotype;

import java.util.Map;
import org.junit.jupiter.api.Test;
import org.pharmgkb.common.util.PathUtils;

import static org.junit.jupiter.api.Assertions.*;


/**
 * JUnit test for {@link VcfReader}.
 *
 * @author Mark Woon
 */
class VcfReaderTest {

  @Test
  void testPhasing() throws Exception {

    VcfReader reader = new VcfReader(PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/VcfReaderTest-phasing.vcf"));
    Map<String, SampleAllele> alleleMap = reader.getAlleleMap();
    assertEquals(6, alleleMap.size());
    for (String chrPos : alleleMap.keySet()) {
      SampleAllele sa = alleleMap.get(chrPos);
      assertEquals(chrPos, sa.getChrPosition());
      if (sa.getChromosome().equals("chr3")) {
        // we treat homogeneous as phased
        assertTrue(sa.isPhased());
      } else {
        assertFalse(sa.isPhased());
      }
    }
  }
}

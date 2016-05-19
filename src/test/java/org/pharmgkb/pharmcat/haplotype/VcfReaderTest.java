package org.pharmgkb.pharmcat.haplotype;

import java.util.Map;
import com.google.common.collect.ImmutableSet;
import org.junit.Test;
import org.pharmgkb.common.util.PathUtils;

import static org.junit.Assert.*;


/**
 * JUnit test for {@link VcfReader}.
 *
 * @author Mark Woon
 */
public class VcfReaderTest {

  @Test
  public void testPhasing() throws Exception {

    ImmutableSet<String> locationsOfInterest = ImmutableSet.of();
    VcfReader reader = new VcfReader(locationsOfInterest,
        PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/VcfReaderTest-phasing.vcf"));
    Map<String, SampleAllele> alleleMap = reader.getAlleleMap();
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

package org.pharmgkb.pharmcat.haplotype;


import java.io.IOException;
import java.nio.file.Path;
import org.junit.jupiter.api.Test;
import org.pharmgkb.common.util.PathUtils;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.hasItems;


/**
 * This is JUnit test for {@link VcfSampleReader}.
 *
 * @author Mark Woon
 */
class VcfSampleReaderTest {

  @Test
  void read() throws IOException {
    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/VcfSampleReaderTest-multisample.vcf");
    VcfSampleReader sampleReader = new VcfSampleReader(vcfFile);

    assertThat(sampleReader.getSamples(), hasItems("Sample_1", "Sample_2"));
  }
}

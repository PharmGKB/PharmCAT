package org.pharmgkb.pharmcat.haplotype;

import java.nio.file.Path;
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

  @Test
  void testFilters() throws Exception {

    DefinitionReader definitionReader = new DefinitionReader();
    Path definitionFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/VcfReaderTest-filters.json");
    definitionReader.read(definitionFile);

    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/VcfReaderTest-filters.vcf");
    VcfReader reader = new VcfReader(NamedAlleleMatcher.calculateLocationsOfInterest(definitionReader), vcfFile);

    assertNotNull(reader.getWarnings());

    assertNotNull(reader.getWarnings().get("chr10:94938828"));
    assertEquals(1, reader.getWarnings().get("chr10:94938828").size());
    assertEquals("Genotype at this position has novel bases (expected G, found C/T in VCF)",
        reader.getWarnings().get("chr10:94938828").iterator().next());

    assertNotNull(reader.getWarnings().get("chr10:94938683"));
    assertEquals(1, reader.getWarnings().get("chr10:94938683").size());
    assertEquals("Discarded genotype at this position because REF in VCF (G) does not match expected reference (A)",
        reader.getWarnings().get("chr10:94938683").iterator().next());

    assertNotNull(reader.getWarnings().get("chr10:94941982"));
    assertEquals(1, reader.getWarnings().get("chr10:94941982").size());
    assertEquals("Discarded genotype at this position because REF in VCF (GA) does not match expected reference (G)",
        reader.getWarnings().get("chr10:94941982").iterator().next());

    assertNotNull(reader.getWarnings().get("chr10:94941897"));
    assertEquals(1, reader.getWarnings().get("chr10:94941897").size());
    assertEquals("PharmCAT preprocessor detected REF mismatch (filter PCATxREF) but this does not match " +
            "current data.  Was the VCF preprocessed with a different version of PharmCAT?",
        reader.getWarnings().get("chr10:94941897").iterator().next());

    assertNotNull(reader.getWarnings().get("chr10:94941915"));
    assertEquals(1, reader.getWarnings().get("chr10:94941915").size());
    assertEquals("PharmCAT preprocessor detected ALT mismatch (filter PCATxALT) but this does not match " +
            "current data.  Was the VCF preprocessed with a different version of PharmCAT?",
        reader.getWarnings().get("chr10:94941915").iterator().next());

    assertEquals(5, reader.getWarnings().size());
  }

  @Test
  void testAdField() throws Exception {

    DefinitionReader definitionReader = new DefinitionReader();
    Path definitionFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/VcfReaderTest-filters.json");
    definitionReader.read(definitionFile);

    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/VcfReaderTest-AD.vcf");
    VcfReader reader = new VcfReader(NamedAlleleMatcher.calculateLocationsOfInterest(definitionReader), vcfFile);

    assertNotNull(reader.getWarnings());

    for (String key : reader.getWarnings().keySet()) {
      System.out.println(key);
      System.out.println("\t" + reader.getWarnings().get(key));
    }

    String pos = "chr10:94942254";
    assertNotNull(reader.getWarnings().get(pos));
    assertEquals(1, reader.getWarnings().get(pos).size());
    assertEquals("Discarding genotype at this position because GT field indicates heterozygous (0/1) but AD field indicates homozygous (91,0)",
        reader.getWarnings().get(pos).iterator().next());

    pos = "chr10:94942255";
    assertNotNull(reader.getWarnings().get(pos));
    assertEquals(1, reader.getWarnings().get(pos).size());
    assertEquals("Discarding genotype at this position because GT field indicates heterozygous (0/1) but AD field indicates homozygous (0,91)",
        reader.getWarnings().get(pos).iterator().next());

    assertEquals(2, reader.getWarnings().size());
  }

  @Test
  void testFiltersNoDefinition() throws Exception {

    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/VcfReaderTest-filters.vcf");
    VcfReader reader = new VcfReader(vcfFile);

    assertNotNull(reader.getWarnings());

    assertNotNull(reader.getWarnings().get("chr10:94938828"));
    assertEquals(1, reader.getWarnings().get("chr10:94938828").size());
    assertTrue(reader.getWarnings().get("chr10:94938828")
        .contains("Genotype at this position has novel bases"));

    assertNotNull(reader.getWarnings().get("chr10:94938683"));
    assertEquals(1, reader.getWarnings().get("chr10:94938683").size());
    assertTrue(reader.getWarnings().get("chr10:94938683")
        .contains("Discarded genotype at this position because REF in VCF (G) does not match expected reference"));

    assertNotNull(reader.getWarnings().get("chr10:94941982"));
    assertEquals(1, reader.getWarnings().get("chr10:94941982").size());
    assertTrue(reader.getWarnings().get("chr10:94941982")
        .contains("Discarded genotype at this position because REF in VCF (GA) does not match expected reference"));

    assertNotNull(reader.getWarnings().get("chr10:94941897"));
    assertEquals(1, reader.getWarnings().get("chr10:94941897").size());
    assertEquals("Discarded genotype at this position because REF in VCF (G) does not match expected reference",
        reader.getWarnings().get("chr10:94941897").iterator().next());

    assertNotNull(reader.getWarnings().get("chr10:94941915"));
    assertEquals(1, reader.getWarnings().get("chr10:94941915").size());
    assertEquals("Genotype at this position has novel bases",
        reader.getWarnings().get("chr10:94941915").iterator().next());

    assertEquals(5, reader.getWarnings().size());
  }
}

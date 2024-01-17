package org.pharmgkb.pharmcat.haplotype;

import java.io.IOException;
import java.nio.file.Path;
import java.util.Map;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.TestInfo;
import org.pharmgkb.common.util.PathUtils;
import org.pharmgkb.parser.vcf.VcfFormatException;
import org.pharmgkb.pharmcat.TestUtils;
import org.pharmgkb.pharmcat.VcfFile;
import org.pharmgkb.pharmcat.definition.DefinitionReader;

import static org.junit.jupiter.api.Assertions.*;


/**
 * JUnit test for {@link VcfReader}.
 *
 * @author Mark Woon
 */
class VcfReaderTest {

  @Test
  void testCompressed() throws Exception {

    VcfReader reader = new VcfReader(PathUtils.getPathToResource("org/pharmgkb/pharmcat/multisample.vcf.bgz"));
    assertEquals(2, reader.getVcfMetadata().getNumSamples());
    Map<String, SampleAllele> alleleMap = reader.getAlleleMap();
    assertEquals(15, alleleMap.size());
  }


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
        assertTrue(sa.isEffectivelyPhased());
      } else {
        assertFalse(sa.isEffectivelyPhased());
      }
    }
  }


  private void printWarnings(VcfReader reader) {
    for (String key : reader.getWarnings().keySet()) {
      System.out.println(key);
      System.out.println("\t" + reader.getWarnings().get(key));
    }
  }


  @Test
  void testGvcfAlt() throws Exception {

    VcfReader reader = new VcfReader(PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/VcfReaderTest-structuralAlt.vcf"));
    assertNotNull(reader.getWarnings());
    printWarnings(reader);

    assertEquals(6, reader.getWarnings().size());
    assertEquals(1, reader.getWarnings().get("chr1:1").size());
    assertFalse((reader.getWarnings().get("chr1:1")).iterator().next().contains("Discarded"));
    assertEquals(1, reader.getWarnings().get("chr1:2").size());
    assertFalse((reader.getWarnings().get("chr1:2")).iterator().next().contains("Discarded"));
    assertEquals(1, reader.getWarnings().get("chr1:3").size());
    assertFalse((reader.getWarnings().get("chr1:3")).iterator().next().contains("Discarded"));
    assertEquals(1, reader.getWarnings().get("chr1:4").size());
    assertTrue((reader.getWarnings().get("chr1:4")).iterator().next().contains("Discarded"));
    assertEquals(1, reader.getWarnings().get("chr1:5").size());
    assertTrue((reader.getWarnings().get("chr1:5")).iterator().next().contains("Discarded"));
    assertEquals(1, reader.getWarnings().get("chr1:6").size());
    assertTrue((reader.getWarnings().get("chr1:6")).iterator().next().contains("Discarded"));

    Map<String, SampleAllele> alleleMap = reader.getAlleleMap();
    assertEquals(3, alleleMap.size());
    assertNotNull(alleleMap.get("chr1:1"));
    assertNotNull(alleleMap.get("chr1:2"));
    assertNotNull(alleleMap.get("chr1:3"));
  }

  @Test
  void testFilters() throws Exception {

    Path definitionFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/VcfReaderTest-filters.json");
    DefinitionReader definitionReader = new DefinitionReader(definitionFile, null);

    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/VcfReaderTest-filters.vcf");
    VcfReader reader = new VcfReader(definitionReader, vcfFile);

    assertNotNull(reader.getWarnings());
    printWarnings(reader);

    assertNotNull(reader.getWarnings().get("chr10:94938828"));
    assertEquals(1, reader.getWarnings().get("chr10:94938828").size());
    assertEquals("The genetic variation at this position does not match what is in the allele definition (expected G, found T in VCF)",
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
            "current data (expected A and got A).  Was the VCF preprocessed with a different version of PharmCAT?",
        reader.getWarnings().get("chr10:94941915").iterator().next());

    assertNotNull(reader.getWarnings().get("chr10:94949281"));
    assertEquals(1, reader.getWarnings().get("chr10:94949281").size());
    assertEquals("ALT uses structural variation '<*>'",
        reader.getWarnings().get("chr10:94949281").iterator().next());

    assertEquals(6, reader.getWarnings().size());
  }

  @Test
  void testAdField() throws Exception {

    Path definitionFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/VcfReaderTest-filters.json");
    DefinitionReader definitionReader = new DefinitionReader(definitionFile, null);

    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/VcfReaderTest-AD.vcf");
    VcfReader reader = new VcfReader(definitionReader, vcfFile);

    assertNotNull(reader.getWarnings());
    printWarnings(reader);

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
  void testAdFieldUnknown() throws Exception {

    Path definitionFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/VcfReaderTest-filters.json");
    DefinitionReader definitionReader = new DefinitionReader(definitionFile, null);

    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/VcfReaderTest-AD-unknown.vcf");
    VcfReader reader = new VcfReader(definitionReader, vcfFile);

    assertNotNull(reader.getWarnings());
    printWarnings(reader);
    assertEquals(0, reader.getWarnings().size());
  }

  @Test
  void testAdFieldInvalid() throws Exception {

    Path definitionFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/VcfReaderTest-filters.json");
    DefinitionReader definitionReader = new DefinitionReader(definitionFile, null);

    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/VcfReaderTest-AD-invalid.vcf");
    VcfReader reader = new VcfReader(definitionReader, vcfFile);

    assertNotNull(reader.getWarnings());
    printWarnings(reader);
    assertEquals(1, reader.getWarnings().size());
    assertEquals(1, reader.getWarnings().get("VCF").size());
    assertEquals("INFO header for AD has unexpected number (f). Expecting 'R'. Treating number as '.' and ignoring AD field.",
        reader.getWarnings().get("VCF").first());
  }

  @Test
  void testAdFieldMissing() throws Exception {

    Path definitionFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/VcfReaderTest-filters.json");
    DefinitionReader definitionReader = new DefinitionReader(definitionFile, null);

    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/VcfReaderTest-AD-missing.vcf");
    VcfReader reader = new VcfReader(definitionReader, vcfFile);

    assertNotNull(reader.getWarnings());
    printWarnings(reader);

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

    assertEquals(3, reader.getWarnings().size());
    assertEquals(1, reader.getWarnings().get("VCF").size());
    assertEquals(VcfReader.MSG_AD_FORMAT_MISSING, reader.getWarnings().get("VCF").first());
  }



  @Test
  void testFiltersNoDefinition() throws Exception {

    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/VcfReaderTest-filters.vcf");
    VcfReader reader = new VcfReader(vcfFile);

    assertNotNull(reader.getWarnings());

    assertNotNull(reader.getWarnings().get("chr10:94938828"));
    assertEquals(1, reader.getWarnings().get("chr10:94938828").size());
    assertTrue(reader.getWarnings().get("chr10:94938828")
        .contains("The genetic variation at this position does not match what is in the allele definition"));

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
    assertEquals("The genetic variation at this position does not match what is in the allele definition",
        reader.getWarnings().get("chr10:94941915").iterator().next());

    assertNotNull(reader.getWarnings().get("chr10:94949281"));
    assertEquals(1, reader.getWarnings().get("chr10:94949281").size());
    assertEquals("ALT uses structural variation '<*>'",
        reader.getWarnings().get("chr10:94949281").iterator().next());

    assertEquals(6, reader.getWarnings().size());
  }


  @Test
  void testIsVcf(TestInfo testInfo) throws IOException {
    Path file = TestUtils.createTempFile(testInfo, ".vcf");
    assertFalse(VcfFile.isVcfFile(file.getParent()));
    assertTrue(VcfFile.isVcfFile(file));
    assertTrue(VcfFile.isVcfFile(TestUtils.createTempFile(testInfo, ".vcf.gz")));
    assertTrue(VcfFile.isVcfFile(TestUtils.createTempFile(testInfo, ".vcf.bgz")));

    assertFalse(VcfFile.isVcfFile(TestUtils.createTempFile(testInfo, ".txt")));
    assertFalse(VcfFile.isVcfFile(TestUtils.createTempFile(testInfo, ".txt.gz")));
  }


  @Test
  void testBadAllele() {

    assertThrows(VcfFormatException.class, () -> {
      VcfReader reader = new VcfReader(PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/VcfReaderTest-badAllele.vcf"));
      assertNotNull(reader.getWarnings());
      printWarnings(reader);
    });
  }
}

package org.pharmgkb.pharmcat.definition;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Optional;
import java.util.Set;
import com.google.common.base.Charsets;
import org.apache.commons.io.FileUtils;
import org.junit.jupiter.api.Test;
import org.pharmgkb.common.util.PathUtils;
import org.pharmgkb.pharmcat.definition.model.DefinitionExemption;
import org.pharmgkb.pharmcat.definition.model.DefinitionFile;
import org.pharmgkb.pharmcat.util.DataSerializer;

import static org.junit.jupiter.api.Assertions.*;


/**
 * JUnit test for {@link DataSerializer}.
 *
 * @author Mark Woon
 */
class DataSerializerTest {

  @Test
  void testJson1() throws Exception {

    // is INS
    Path inFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/definition/CYP3A5.good.tsv");
    DefinitionFile[] definitionFiles = testJson(inFile);
    CuratedDefinitionParserTest.assertInsertFromCyp3a5GoodTsv(definitionFiles[1]);
  }


  @Test
  void testJsonRepeat() throws Exception {

    Path inFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/definition/repeats.tsv");
    DefinitionFile[] definitionFiles = testJson(inFile);
    CuratedDefinitionParserTest.assertRepeatFromRepeatsTsv(definitionFiles[1]);
  }


  @Test
  void testJson2() throws Exception {

    // contains population frequencies
    Path inFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/definition/CYP2C19.tsv");
    DefinitionFile[] definitionFiles = testJson(inFile);
    assertEquals(8, definitionFiles[1].getPopulations().size());
    assertNotNull(definitionFiles[1].getNamedAlleles().get(0).getPopFreqMap());
    assertFalse(definitionFiles[1].getNamedAlleles().get(0).getPopFreqMap().isEmpty());
    assertEquals("0.331", definitionFiles[1].getNamedAlleles().get(0).getPopFreqMap().get("African Allele Frequency"));
  }


  private DefinitionFile[] testJson(Path inFile) throws IOException {
    DefinitionFile definitionFile = new CuratedDefinitionParser(inFile).parse();
    assertNotNull(definitionFile.getVariants());
    assertNotNull(definitionFile.getVariantAlleles());
    assertEquals(definitionFile.getVariants().length, definitionFile.getVariantAlleles().size());

    Path jsonFile = inFile.getParent().resolve(PathUtils.getBaseFilename(inFile) + ".json");
    DataSerializer cdSerializer = new DataSerializer();
    // write it out
    cdSerializer.serializeToJson(definitionFile, jsonFile);
    // read it back in
    DefinitionFile jsonDefinitionFile = cdSerializer.deserializeDefinitionsFromJson(jsonFile);
    assertNotNull(jsonDefinitionFile.getVariants());
    assertNotNull(jsonDefinitionFile.getVariantAlleles());
    assertEquals(jsonDefinitionFile.getVariants().length, jsonDefinitionFile.getVariantAlleles().size());

    assertEquals(definitionFile, jsonDefinitionFile);
    return new DefinitionFile[] { definitionFile, jsonDefinitionFile };
  }


  @Test
  void testExemptions() throws Exception {

    Path tsvFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/definition/exemptions.tsv");
    Path jsonFile = Files.createTempFile("transformExemptions", ".json");
    Path refJsonFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/definition/exemptions.json");

    try {
      DataSerializer dataSerializer = new DataSerializer();
      Set<DefinitionExemption> definitionExemptions = dataSerializer.deserializeExemptionsFromTsv(tsvFile);
      assertEquals(3, definitionExemptions.size());

      Optional<DefinitionExemption> cyp2c9 = definitionExemptions.stream()
          .filter((de) -> de.getGene().equals("CYP2C9"))
          .findFirst();
      assertTrue(cyp2c9.isPresent());
      assertEquals(0, cyp2c9.get().getIgnoredPositions().size());
      assertEquals(1, cyp2c9.get().getExtraPositions().size());
      assertEquals(0, cyp2c9.get().getIgnoredAlleles().size());
      assertFalse(cyp2c9.get().isAllHits());
      assertTrue(cyp2c9.get().isAssumeReference());

      Optional<DefinitionExemption> cyp2c19 = definitionExemptions.stream()
          .filter((de) -> de.getGene().equals("CYP2C19"))
          .findFirst();
      assertTrue(cyp2c19.isPresent());
      assertEquals(1, cyp2c19.get().getIgnoredPositions().size());
      assertEquals(0, cyp2c19.get().getExtraPositions().size());
      assertEquals(0, cyp2c19.get().getIgnoredAlleles().size());
      assertFalse(cyp2c19.get().isAllHits());
      assertTrue(cyp2c19.get().isAssumeReference());

      Optional<DefinitionExemption> ugt1a1 = definitionExemptions.stream()
          .filter((de) -> de.getGene().equals("UGT1A1"))
          .findFirst();
      assertTrue(ugt1a1.isPresent());
      assertEquals(0, ugt1a1.get().getIgnoredPositions().size());
      assertEquals(0, ugt1a1.get().getExtraPositions().size());
      assertEquals(0, ugt1a1.get().getIgnoredAlleles().size());
      assertTrue(ugt1a1.get().isAllHits());
      assertFalse(ugt1a1.get().isAssumeReference());

      // write it out
      dataSerializer.serializeToJson(definitionExemptions, jsonFile);
      // compare with expected
      FileUtils.contentEqualsIgnoreEOL(refJsonFile.toFile(), jsonFile.toFile(), Charsets.UTF_8.displayName());


    } finally {
      Files.deleteIfExists(jsonFile);
    }
  }
}

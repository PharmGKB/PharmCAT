package org.pharmgkb.pharmcat.definition;

import java.io.IOException;
import java.nio.file.Path;
import org.junit.jupiter.api.Test;
import org.pharmgkb.common.util.PathUtils;
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
}

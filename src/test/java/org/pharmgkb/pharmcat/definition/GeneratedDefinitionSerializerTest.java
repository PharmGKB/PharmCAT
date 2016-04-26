package org.pharmgkb.pharmcat.definition;

import java.io.IOException;
import java.nio.file.Path;
import org.junit.Test;
import org.pharmgkb.common.util.PathUtils;
import org.pharmgkb.pharmcat.TestUtil;
import org.pharmgkb.pharmcat.definition.model.DefinitionFile;

import static org.junit.Assert.assertEquals;


/**
 * JUnit test for {@link GeneratedDefinitionSerializer}.
 *
 * @author Mark Woon
 */
public class GeneratedDefinitionSerializerTest {

  @Test
  public void testJson1() throws Exception {

    Path inFile = TestUtil.getFile("org/pharmgkb/pharmcat/definition/CYP3A5.good.tsv");
    testJson(inFile);
  }

  @Test
  public void testJson2() throws Exception {

    // contains population frequencies
    Path inFile = TestUtil.getFile("org/pharmgkb/pharmcat/definition/CYP2C19.tsv");
    testJson(inFile);
  }

  private void testJson(Path inFile) throws IOException {
    DefinitionFile definitionFile = new CuratedDefinitionParser(inFile).parse();

    Path jsonFile = inFile.getParent().resolve(PathUtils.getBaseFilename(inFile) + ".json");
    GeneratedDefinitionSerializer cdSerializer = new GeneratedDefinitionSerializer();
    // write it out
    cdSerializer.serializeToJson(definitionFile, jsonFile);
    // read it back in
    DefinitionFile jsonDefinitionFile = cdSerializer.deserializeFromJson(jsonFile);

    assertEquals(definitionFile, jsonDefinitionFile);
  }


  @Test
  public void testTsv() throws Exception {

    Path inFile = TestUtil.getFile("org/pharmgkb/pharmcat/definition/CYP3A5.good.tsv");
    DefinitionFile definitionFile = new CuratedDefinitionParser(inFile).parse();

    Path tsvFile = inFile.getParent().resolve(PathUtils.getBaseFilename(inFile) + "-generated.tsv");
    GeneratedDefinitionSerializer cdSerializer = new GeneratedDefinitionSerializer();
    // write it out
    cdSerializer.serializeToTsv(definitionFile, tsvFile);

    // TODO(markwoon): check that file is correct!
    System.out.println(tsvFile);
  }
}

package org.pharmgkb.pharmcat.definition;

import java.nio.file.Path;
import com.google.common.base.Joiner;
import org.junit.jupiter.api.Test;
import org.pharmgkb.common.util.PathUtils;
import org.pharmgkb.pharmcat.ParseException;
import org.pharmgkb.pharmcat.definition.model.DefinitionFile;
import org.pharmgkb.pharmcat.definition.model.NamedAllele;
import org.pharmgkb.pharmcat.definition.model.VariantType;

import static org.junit.jupiter.api.Assertions.*;


class CuratedDefinitionParserTest {

  @Test
  void testReaderGood() {

    Path tsvFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/definition/CYP3A5.good.tsv");
    CuratedDefinitionParser parser = new CuratedDefinitionParser(tsvFile);
    assertInsertFromCyp3a5GoodTsv(parser.parse());
  }

  static void assertInsertFromCyp3a5GoodTsv(DefinitionFile definitionFile) {
    assertEquals(9, definitionFile.getNamedAlleles().size());
    assertEquals(8, definitionFile.getVariants().length);
    assertEquals("g.99652770_99652771insA", definitionFile.getVariants()[6].getChromosomeHgvsName());
    assertEquals(99652770, definitionFile.getVariants()[6].getPosition());
    assertEquals(99652770, definitionFile.getVariants()[6].getVcfPosition());
    assertEquals(VariantType.INS, definitionFile.getVariants()[6].getType());
  }


  @Test
  void testReader2C19() {

    Path tsvFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/definition/CYP2C19.tsv");
    CuratedDefinitionParser parser = new CuratedDefinitionParser(tsvFile);
    DefinitionFile definitionFile = parser.parse();
    assertEquals("CYP2C19", definitionFile.getGeneSymbol());
    assertEquals(38, definitionFile.getVariants().length);
  }

  @Test
  void testReaderCFTR() {
    
    Path tsvFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/definition/CFTR.tsv");
    CuratedDefinitionParser parser = new CuratedDefinitionParser(tsvFile);
    DefinitionFile definitionFile = parser.parse();
    assertEquals("CFTR", definitionFile.getGeneSymbol());
    assertEquals(34, definitionFile.getVariants().length);
    assertEquals(10, definitionFile.getNotes().size());
  }

  @Test
  void testReaderRepeats() {

    Path tsvFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/definition/repeats.tsv");
    CuratedDefinitionParser parser = new CuratedDefinitionParser(tsvFile);
    assertRepeatFromRepeatsTsv(parser.parse());
  }

  static void assertRepeatFromRepeatsTsv(DefinitionFile definitionFile) {
    assertEquals(8, definitionFile.getNamedAlleles().size());
    assertEquals(233760234, definitionFile.getVariants()[2].getPosition());
    assertEquals(233760234, definitionFile.getVariants()[2].getVcfPosition());
    assertEquals(VariantType.REPEAT, definitionFile.getVariants()[2].getType());
    assertEquals("A(TA)6TAA", definitionFile.getVariants()[2].getReferenceRepeat());

    NamedAllele na = definitionFile.getNamedAlleles().iterator().next();
    System.out.println(Joiner.on("; ").join(na.getAlleles()));
  }

  @Test
  void testReaderDeletes() {
    // TODO(markwoon): finish this
  }



  @Test
  void testReaderBad() {

    try {
      Path tsvFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/definition/CYP3A5.bad.tsv");
      CuratedDefinitionParser parser = new CuratedDefinitionParser(tsvFile);
      DefinitionFile definitionFile = parser.parse();
      assertEquals(9, definitionFile.getNamedAlleles().size());
      fail("Missed invalid alleles");
    } catch (ParseException ex) {
      System.out.println(ex.getMessage());
      assertTrue(ex.getMessage().contains("Invalid bases (BBQ)"));
    }
  }

  @Test
  void testColumnNumberToName() {
    assertEquals("A", CuratedDefinitionParser.columnNumberToName(0));
    assertEquals("B", CuratedDefinitionParser.columnNumberToName(1));
    assertEquals("Z", CuratedDefinitionParser.columnNumberToName(25));
    assertEquals("AA", CuratedDefinitionParser.columnNumberToName(26));
    assertEquals("AB", CuratedDefinitionParser.columnNumberToName(27));
    assertEquals("AZ", CuratedDefinitionParser.columnNumberToName(51));
    assertEquals("BA", CuratedDefinitionParser.columnNumberToName(52));
  }
}

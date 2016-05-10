package org.pharmgkb.pharmcat.definition;

import java.nio.file.Path;
import com.google.common.base.Joiner;
import org.junit.Test;
import org.pharmgkb.pharmcat.ParseException;
import org.pharmgkb.pharmcat.TestUtil;
import org.pharmgkb.pharmcat.definition.model.DefinitionFile;
import org.pharmgkb.pharmcat.definition.model.NamedAllele;
import org.pharmgkb.pharmcat.definition.model.VariantType;

import static org.junit.Assert.*;


public class CuratedDefinitionParserTest {

  @Test
  public void testReaderGood() throws Exception {

    Path tsvFile = TestUtil.getFile("org/pharmgkb/pharmcat/definition/CYP3A5.good.tsv");
    CuratedDefinitionParser parser = new CuratedDefinitionParser(tsvFile);
    DefinitionFile definitionFile = parser.parse();
    assertEquals(9, definitionFile.getNamedAlleles().size());
    NamedAllele na = definitionFile.getNamedAlleles().iterator().next();
    System.out.println(Joiner.on("; ").join(na.getAlleles()));
  }


  @Test
  public void testReaderRepeats() throws Exception {

    Path tsvFile = TestUtil.getFile("org/pharmgkb/pharmcat/definition/repeats.tsv");
    CuratedDefinitionParser parser = new CuratedDefinitionParser(tsvFile);
    DefinitionFile definitionFile = parser.parse();
    assertEquals(8, definitionFile.getNamedAlleles().size());
    NamedAllele na = definitionFile.getNamedAlleles().iterator().next();
    assertEquals(233760234, definitionFile.getVariants()[2].getPosition());
    assertEquals(VariantType.REPEAT, definitionFile.getVariants()[2].getType());
    assertEquals("A(TA)6TAA", definitionFile.getVariants()[2].getReferenceRepeat());
    System.out.println(Joiner.on("; ").join(na.getAlleles()));
  }


  @Test
  public void testReaderBad() throws Exception {

    try {
      Path tsvFile = TestUtil.getFile("org/pharmgkb/pharmcat/definition/CYP3A5.bad.tsv");
      CuratedDefinitionParser parser = new CuratedDefinitionParser(tsvFile);
      DefinitionFile definitionFile = parser.parse();
      assertEquals(9, definitionFile.getNamedAlleles().size());
      fail("Missed invalid alleles");
    } catch (ParseException ex) {
      System.out.println(ex.getMessage());
      assertTrue(ex.getMessage().contains("Invalid bases (BBQ)"));
    }
  }
}

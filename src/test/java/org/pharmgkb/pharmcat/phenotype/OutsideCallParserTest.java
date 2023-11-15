package org.pharmgkb.pharmcat.phenotype;

import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Path;
import java.util.List;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.TestInfo;
import org.pharmgkb.pharmcat.TestUtils;
import org.pharmgkb.pharmcat.phenotype.model.OutsideCall;
import org.pharmgkb.pharmcat.reporter.BadOutsideCallException;

import static org.junit.jupiter.api.Assertions.*;


class OutsideCallParserTest {

  @Test
  void testMinimalInput(TestInfo testInfo) throws IOException {
    Path outsideCallPath = TestUtils.createTestFile(testInfo, ".tsv");
    try (FileWriter fw = new FileWriter(outsideCallPath.toFile())) {
      fw.write("CYP2C9\t*1/*2");
    }

    List<OutsideCall> calls = OutsideCallParser.parse(outsideCallPath);
    assertNotNull(calls);
    assertEquals(1, calls.size());
    assertEquals("CYP2C9", calls.get(0).getGene());
    assertEquals("*1/*2", calls.get(0).getDiplotype());
  }

  @Test
  void testTwoGenes(TestInfo testInfo) throws IOException {
    Path tmpOutsideCallPath = TestUtils.createTestFile(testInfo, ".tsv");
    try (FileWriter fw = new FileWriter(tmpOutsideCallPath.toFile())) {
      fw.write("""
          CYP2C9\t*1/*2
          CYP2C19\t*3/*4
          """);
    }

    List<OutsideCall> calls = OutsideCallParser.parse(tmpOutsideCallPath);
    assertNotNull(calls);
    assertEquals(2, calls.size());
    assertEquals("CYP2C9", calls.get(0).getGene());
    assertEquals("*1/*2", calls.get(0).getDiplotype());

    assertEquals("CYP2C19", calls.get(1).getGene());
    assertEquals("*3/*4", calls.get(1).getDiplotype());
  }

  @Test
  void testBadFormat(TestInfo testInfo) throws IOException {
    Path outsideCallPath = TestUtils.createTestFile(testInfo, ".tsv");
    try (FileWriter fw = new FileWriter(outsideCallPath.toFile())) {
      fw.write("""
          CYP2C9\t*1/*2
          CYP2C19\t*3/*4/*2
          """);
    }

    assertThrows(BadOutsideCallException.class, () -> OutsideCallParser.parse(outsideCallPath));
  }


  @Test
  void testCommentsAndEmptyLines(TestInfo testInfo) throws IOException {
    Path outsideCallPath = TestUtils.createTestFile(testInfo, ".tsv");
    try (FileWriter fw = new FileWriter(outsideCallPath.toFile())) {
      fw.write("""
          # comment
          CYP2C9\t*1/*2
          
          ## another comment
          
          CYP2C9\t*3/*4
          
          """);
    }

    List<OutsideCall> calls = OutsideCallParser.parse(outsideCallPath);
    assertEquals(2, calls.size());
    assertEquals("*1/*2", calls.get(0).getDiplotype());
    assertEquals("*3/*4", calls.get(1).getDiplotype());
  }
}

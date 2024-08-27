package org.pharmgkb.pharmcat.phenotype;

import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Path;
import java.util.List;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.TestInfo;
import org.pharmgkb.pharmcat.Env;
import org.pharmgkb.pharmcat.TestUtils;
import org.pharmgkb.pharmcat.phenotype.model.OutsideCall;
import org.pharmgkb.pharmcat.reporter.BadOutsideCallException;

import static org.junit.jupiter.api.Assertions.*;


class OutsideCallParserTest {
  private static Env s_env;

  @BeforeAll
  static void prepare() throws Exception {
    s_env = new Env();
  }


  @Test
  void testMinimalInput(TestInfo testInfo) throws IOException {
    Path outsideCallPath = TestUtils.createTestFile(testInfo, ".tsv");
    try (FileWriter fw = new FileWriter(outsideCallPath.toFile())) {
      fw.write("CYP2C9\t*1/*2");
    }

    List<OutsideCall> calls = OutsideCallParser.parse(s_env, outsideCallPath);
    assertNotNull(calls);
    assertEquals(1, calls.size());
    assertEquals("CYP2C9", calls.get(0).getGene());
    assertEquals("*1/*2", calls.get(0).getDiplotype());
  }

  @Test
  void testGenePrefixStripping(TestInfo testInfo) throws IOException {
    Path outsideCallPath = TestUtils.createTestFile(testInfo, ".tsv");
    try (FileWriter fw = new FileWriter(outsideCallPath.toFile())) {
      fw.write("CYP2C9\tCYP2C9*1/CYP2C9      *2");
    }

    List<OutsideCall> calls = OutsideCallParser.parse(s_env, outsideCallPath);
    assertNotNull(calls);
    assertEquals(1, calls.size());
    assertEquals("CYP2C9", calls.get(0).getGene());
    assertEquals("*1/*2", calls.get(0).getDiplotype());
  }

  @Test
  void testPhenotype(TestInfo testInfo) throws IOException {
    Path outsideCallPath = TestUtils.createTestFile(testInfo, ".tsv");
    try (FileWriter fw = new FileWriter(outsideCallPath.toFile())) {
      fw.write("CYP2C9\t*1/*2\tNormal Metabolizer");
    }

    List<OutsideCall> calls = OutsideCallParser.parse(s_env, outsideCallPath);
    assertNotNull(calls);
    assertEquals(1, calls.size());
    assertEquals("CYP2C9", calls.get(0).getGene());
    assertEquals("*1/*2", calls.get(0).getDiplotype());
    assertEquals("Normal Metabolizer", calls.get(0).getPhenotype());
  }

  @Test
  void testPrefixedPhenotype(TestInfo testInfo) throws IOException {
    Path outsideCallPath = TestUtils.createTestFile(testInfo, ".tsv");
    try (FileWriter fw = new FileWriter(outsideCallPath.toFile())) {
      fw.write("CYP2C9\t*1/*2\tCYP2C9 Normal Metabolizer");
    }

    List<OutsideCall> calls = OutsideCallParser.parse(s_env, outsideCallPath);
    assertNotNull(calls);
    assertEquals(1, calls.size());
    assertEquals("CYP2C9", calls.get(0).getGene());
    assertEquals("*1/*2", calls.get(0).getDiplotype());
    assertEquals("Normal Metabolizer", calls.get(0).getPhenotype());
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

    List<OutsideCall> calls = OutsideCallParser.parse(s_env, tmpOutsideCallPath);
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

    assertThrows(BadOutsideCallException.class, () -> OutsideCallParser.parse(s_env, outsideCallPath));
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

    List<OutsideCall> calls = OutsideCallParser.parse(s_env, outsideCallPath);
    assertEquals(2, calls.size());
    assertEquals("*1/*2", calls.get(0).getDiplotype());
    assertEquals("*3/*4", calls.get(1).getDiplotype());
  }
}

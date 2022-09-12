package org.pharmgkb.pharmcat.reporter.io;

import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Path;
import java.util.List;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.TestInfo;
import org.pharmgkb.pharmcat.ParseException;
import org.pharmgkb.pharmcat.TestUtils;
import org.pharmgkb.pharmcat.reporter.model.OutsideCall;

import static org.junit.jupiter.api.Assertions.*;


class OutsideCallParserTest {

  @Test
  void testAstrolabeFormat(TestInfo testInfo) throws IOException {
    String sf_astrolabeOutput = "##Test Astrolabe output\n" +
        "#ROI_label\tdiplotype labels\tdiplotype activity\tdiplotype calling notes\tjaccard\tpart\tpValue\tROI notes\t" +
        "special case\tnomenclature version\n" +
        "CYP2D6\tCYP2D6*1/CYP2D6*4\t\t\t0.6\t0.75\tp: 0.0\t\t\tv1.9-2017_02_09\n";

    Path tempAstroPath = TestUtils.createTestFile(testInfo, ".tsv");
    try (FileWriter fw = new FileWriter(tempAstroPath.toFile())) {
      fw.write(sf_astrolabeOutput);
    }

    List<OutsideCall> calls = OutsideCallParser.parse(tempAstroPath);
    assertNotNull(calls);
    assertEquals(1, calls.size());
    assertEquals("CYP2D6", calls.get(0).getGene());
    assertEquals(1, calls.get(0).getDiplotypes().size());
    assertEquals("*1/*4", calls.get(0).getDiplotypes().get(0));
  }

  @Test
  void testAstrolabeMultiple(TestInfo testInfo) throws IOException {
    String sf_astrolabeOutput = "##Test Astrolabe output\n" +
        "#ROI_label\tdiplotype labels\tdiplotype activity\tdiplotype calling notes\tjaccard\tpart\tpValue\tROI notes\t" +
        "special case\tnomenclature version\n" +
        "CYP2D6\tCYP2D6*1/CYP2D6*4 or *2/*3\t\t\t0.6\t0.75\tp: 0.0\t\t\tv1.9-2017_02_09\n";

    Path tempAstroPath = TestUtils.createTestFile(testInfo, ".tsv");
    try (FileWriter fw = new FileWriter(tempAstroPath.toFile())) {
      fw.write(sf_astrolabeOutput);
    }

    List<OutsideCall> calls = OutsideCallParser.parse(tempAstroPath);
    assertNotNull(calls);
    assertEquals(1, calls.size());
    assertEquals("CYP2D6", calls.get(0).getGene());
    assertEquals(2, calls.get(0).getDiplotypes().size());
    assertTrue(calls.get(0).getDiplotypes().contains("*1/*4"));
    assertTrue(calls.get(0).getDiplotypes().contains("*2/*3"));
  }

  @Test
  void testMinimalInput(TestInfo testInfo) throws IOException {
    String sf_astrolabeOutput = "CYP2C9\t*1/*2";

    Path tempAstroPath = TestUtils.createTestFile(testInfo, ".tsv");
    try (FileWriter fw = new FileWriter(tempAstroPath.toFile())) {
      fw.write(sf_astrolabeOutput);
    }

    List<OutsideCall> calls = OutsideCallParser.parse(tempAstroPath);
    assertNotNull(calls);
    assertEquals(1, calls.size());
    assertEquals("CYP2C9", calls.get(0).getGene());
    assertEquals(1, calls.get(0).getDiplotypes().size());
    assertEquals("*1/*2", calls.get(0).getDiplotypes().get(0));
  }

  @Test
  void testTwoGenes(TestInfo testInfo) throws IOException {
    String sf_astrolabeOutput = "CYP2C9\t*1/*2\nCYP2C19\t*3/*4";

    Path tempAstroPath = TestUtils.createTestFile(testInfo, ".tsv");
    try (FileWriter fw = new FileWriter(tempAstroPath.toFile())) {
      fw.write(sf_astrolabeOutput);
    }

    List<OutsideCall> calls = OutsideCallParser.parse(tempAstroPath);
    assertNotNull(calls);
    assertEquals(2, calls.size());
    assertEquals("CYP2C9", calls.get(0).getGene());
    assertEquals(1, calls.get(0).getDiplotypes().size());
    assertEquals("*1/*2", calls.get(0).getDiplotypes().get(0));

    assertEquals("CYP2C19", calls.get(1).getGene());
    assertEquals(1, calls.get(1).getDiplotypes().size());
    assertEquals("*3/*4", calls.get(1).getDiplotypes().get(0));
  }

  @Test
  void testBadFormat(TestInfo testInfo) throws IOException {
    String sf_astrolabeOutput = "CYP2C9\t*1/*2\nCYP2C19\t*3/*4/*2";

    Path tempAstroPath = TestUtils.createTestFile(testInfo, ".tsv");
    try (FileWriter fw = new FileWriter(tempAstroPath.toFile())) {
      fw.write(sf_astrolabeOutput);
    }

    try {
      OutsideCallParser.parse(tempAstroPath);
      fail("Bad ouside call format should have failed");
    }
    catch (ParseException ex) {
      // ignore this exception since we want it to happen
    }
  }
}

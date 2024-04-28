package org.pharmgkb.pharmcat.haplotype;

import java.io.IOException;
import java.nio.file.Path;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.TestInfo;
import org.pharmgkb.pharmcat.TestUtils;
import org.pharmgkb.pharmcat.TestVcfBuilder;
import org.pharmgkb.pharmcat.definition.DefinitionReader;
import org.pharmgkb.pharmcat.definition.model.DefinitionFile;
import org.pharmgkb.pharmcat.haplotype.model.DiplotypeMatch;
import org.pharmgkb.pharmcat.haplotype.model.GeneCall;
import org.pharmgkb.pharmcat.haplotype.model.Result;
import org.pharmgkb.pharmcat.util.DataManager;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.pharmgkb.pharmcat.haplotype.NamedAlleleMatcherTest.testMatchNamedAlleles;


/**
 * RYR1-specific tests for {@link NamedAlleleMatcher}.
 *
 * @author Mark Woon
 */
public class NamedAlleleMatcherRyr1Test {
  private static final Path sf_definitionFile = DataManager.DEFAULT_DEFINITION_DIR.resolve("RYR1_translation.json");
  private static DefinitionFile s_definitionFile;


  @BeforeAll
  static void prepare() throws IOException {
    //TestUtils.setSaveTestOutput(true);
    DefinitionReader definitionReader = new DefinitionReader(sf_definitionFile, null);
    s_definitionFile = definitionReader.getDefinitionFile("RYR1");
  }

  @AfterEach
  void deleteDirectory(TestInfo testInfo) {
    TestUtils.deleteTestOutputDirectory(testInfo);
  }



  @Test
  void reference(TestInfo testInfo) throws Exception {

    Result rez = testMatchNamedAlleles(sf_definitionFile, new TestVcfBuilder(testInfo, "ref/ref")
        .reference("RYR1")
        .generate());
    assertEquals(1, rez.getGeneCalls().size());
    GeneCall gc = rez.getGeneCalls().get(0);
    assertEquals("RYR1", gc.getGene());
    assertEquals(1, gc.getDiplotypes().size());
    DiplotypeMatch dm = gc.getDiplotypes().iterator().next();
    assertEquals("Reference/Reference", dm.toString());

    // testing to make sure variants and sequences get set correctly
    assertEquals(s_definitionFile.getVariants().length, gc.getVariants().size());

    assertEquals(1, dm.getHaplotype1().getSequences().size());
    assertNotNull(dm.getHaplotype2());
    assertEquals(1, dm.getHaplotype2().getSequences().size());
    assertEquals(dm.getHaplotype1().getSequences().first(), dm.getHaplotype2().getSequences().first());
    assertEquals(1, dm.getSequences().size());
    String[] seqPair = dm.getSequences().iterator().next();
    assertEquals(seqPair[0], seqPair[1]);
    assertEquals(dm.getHaplotype1().getSequences().first(), seqPair[0]);
  }


  @Test
  void rs137933390Het(TestInfo testInfo) throws Exception {

    Result rez = testMatchNamedAlleles(sf_definitionFile, new TestVcfBuilder(testInfo, "ref/ref")
        .variation("RYR1", "rs137933390", "A", "G")  // c.4178A>G
        .generate());

    assertEquals(1, rez.getGeneCalls().size());
    GeneCall gc = rez.getGeneCalls().get(0);
    assertEquals("RYR1", gc.getGene());
    assertEquals(1, gc.getDiplotypes().size());
    DiplotypeMatch dm = gc.getDiplotypes().iterator().next();
    assertEquals("Reference/c.4178A>G", dm.toString());

    // testing to make sure variants and sequences get set correctly
    assertEquals(s_definitionFile.getVariants().length, gc.getVariants().size());
    assertEquals(0, gc.getMatchData().getMissingPositions().size());

    assertEquals(1, dm.getHaplotype1().getSequences().size());
    assertNotNull(dm.getHaplotype2());
    assertEquals(1, dm.getHaplotype2().getSequences().size());
    assertEquals(1, dm.getSequences().size());
  }

  @Test
  void rs137933390Het_missing1(TestInfo testInfo) throws Exception {

    Result rez = testMatchNamedAlleles(sf_definitionFile, new TestVcfBuilder(testInfo, "ref/ref")
        .variation("RYR1", "rs137933390", "A", "G")  // c.4178A>G
        .missing("RYR1", "rs193922753")
        .generate());

    assertEquals(1, rez.getGeneCalls().size());
    GeneCall gc = rez.getGeneCalls().get(0);
    assertEquals("RYR1", gc.getGene());
    assertEquals(1, gc.getDiplotypes().size());
    DiplotypeMatch dm = gc.getDiplotypes().iterator().next();
    assertEquals("Reference/c.4178A>G", dm.toString());

    // testing to make sure variants and sequences get set correctly
    assertEquals(s_definitionFile.getVariants().length, gc.getVariants().size() + 1);
    assertEquals(1, gc.getMatchData().getMissingPositions().size());

    assertEquals(1, dm.getHaplotype1().getSequences().size());
    assertNotNull(dm.getHaplotype2());
    assertEquals(1, dm.getHaplotype2().getSequences().size());
    assertEquals(1, dm.getSequences().size());
  }
}

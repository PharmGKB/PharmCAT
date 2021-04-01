package org.pharmgkb.pharmcat.phenotype;

import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.function.Supplier;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import org.pharmgkb.pharmcat.VcfTestUtils;
import org.pharmgkb.pharmcat.haplotype.DefinitionReader;
import org.pharmgkb.pharmcat.haplotype.NamedAlleleMatcher;
import org.pharmgkb.pharmcat.haplotype.model.GeneCall;
import org.pharmgkb.pharmcat.haplotype.model.Result;
import org.pharmgkb.pharmcat.reporter.io.OutsideCallParser;
import org.pharmgkb.pharmcat.reporter.model.result.GeneReport;
import org.pharmgkb.pharmcat.util.DataManager;

import static org.junit.jupiter.api.Assertions.*;


/**
 * Unit test class for the {@link Phenotyper} class. Will run through synthetic data and assert genes are found and
 * diplotypes are called with certain results.
 */
class PhenotyperTest {
  private static final Supplier<RuntimeException> unfoundGene = () -> new RuntimeException("Gene report not found");
  private static DefinitionReader s_definitionReader;

  @BeforeAll
  static void setup() throws IOException {
    s_definitionReader = new DefinitionReader();
    s_definitionReader.read(DataManager.DEFAULT_DEFINITION_DIR);
  }

  @Test
  void testCyp2C19Het() throws Exception {
    Phenotyper phenotyper = new Phenotyper(
        matchVcfData("cyp2c19/s4s17het.vcf"),
        OutsideCallParser.parse("CYP2D6\t*1/*3"));

    assertCalledByMatcher(phenotyper, "CYP2C19", "CYP2D6");
    assertNotCalledByMatcher(phenotyper, "CYP2C9");

    assertDiplotypeDisplay(phenotyper, "CYP2C19", "*1/*4");
    assertLookup(phenotyper, "CYP2C19", "*1", "*4");

    assertDiplotypeDisplay(phenotyper, "CYP2D6", "*1/*3");
    assertLookup(phenotyper, "CYP2D6", "*1", "*3");
  }

  @Test
  void testCyp2D6Only() {
    Phenotyper phenotyper = new Phenotyper(
        new ArrayList<>(),
        OutsideCallParser.parse("CYP2D6\t*1/*3"));

    assertCalledByMatcher(phenotyper, "CYP2D6");

    assertFalse(phenotyper.findGeneReport("CYP2C9").isPresent(),
        "CYP2C9 report should not be present since the NamedAlleleMatcher was never called");

    assertDiplotypeDisplay(phenotyper, "CYP2D6", "*1/*3");
  }

  @Test
  void testCyp2C19Hom() throws Exception {
    Phenotyper phenotyper = new Phenotyper(
        matchVcfData("cyp2c19/s2s2.vcf"),
        new ArrayList<>());

    assertCalledByMatcher(phenotyper, "CYP2C19");
    assertNotCalledByMatcher(phenotyper, "CYP2D6");

    assertDiplotypeDisplay(phenotyper, "CYP2C19", "*2/*2");
  }

  @Test
  void testUGT1A1Phased() throws Exception {
    Phenotyper phenotyper = new Phenotyper(
        matchVcfData("UGT1A1/s1s60s80phased.vcf"),
        new ArrayList<>());

    assertCalledByMatcher(phenotyper, "UGT1A1");

    assertDiplotypeDisplay(phenotyper, "UGT1A1", "*1/*80");
    assertLookup(phenotyper, "UGT1A1", "*1", "*80");
    assertTrue(phenotyper.findGeneReport("UGT1A1").orElseThrow(unfoundGene).isPhased());
  }

  @Test
  void testUGT1A1Unphased() throws Exception {
    Phenotyper phenotyper = new Phenotyper(
        matchVcfData("UGT1A1/s1s60s80phased.vcf"),
        new ArrayList<>());

    assertCalledByMatcher(phenotyper, "UGT1A1");

    assertDiplotypeDisplay(phenotyper, "UGT1A1", "*1/*80");
    assertLookup(phenotyper, "UGT1A1", "*1", "*80");
    assertTrue(phenotyper.findGeneReport("UGT1A1").orElseThrow(unfoundGene).isPhased());
  }

  @Test
  void testNUDT15() throws Exception {
    Phenotyper phenotyper = new Phenotyper(
        matchVcfData("NUDT15/refref.vcf"),
        new ArrayList<>());

    assertCalledByMatcher(phenotyper, "NUDT15");

    assertDiplotypeDisplay(phenotyper, "NUDT15", "*1/*1");
    assertLookup(phenotyper, "NUDT15", "*1", "*1");
  }

  @Test
  void testNUDT15star3() throws Exception {
    Phenotyper phenotyper = new Phenotyper(
        matchVcfData("NUDT15/s3ref.vcf"),
        new ArrayList<>());

    assertCalledByMatcher(phenotyper, "NUDT15");

    assertDiplotypeDisplay(phenotyper, "NUDT15", "*1/*3");
    assertLookup(phenotyper, "NUDT15", "*1", "*3");
  }


  //  Helper methods found below =======================================================================================

  /**
   * Run the {@link NamedAlleleMatcher} on the specified VCF data and return the {@link GeneCall} results
   * @param vcfs VCF test fils found in the test resources package
   * @return a List of GeneCall objects
   */
  private List<GeneCall> matchVcfData(String ...vcfs) throws IOException {
    NamedAlleleMatcher matcher = new NamedAlleleMatcher(s_definitionReader);

    Path tempVcfPath = Files.createTempFile(getClass().getSimpleName(), ".vcf");
    try (FileWriter fw = new FileWriter(tempVcfPath.toFile())) {
      fw.write(VcfTestUtils.writeVcf(vcfs));
    } catch (Exception ex) {
      ex.printStackTrace();
      throw ex;
    }

    VcfTestUtils.writeVcf(vcfs);
    Result result = matcher.call(tempVcfPath);
    return result.getGeneCalls();
  }

  /**
   * Check to see if all the given genes have been called by the matcher
   */
  private void assertCalledByMatcher(Phenotyper phenotyper, String... genes) {
    Arrays.stream(genes)
        .forEach(g -> assertTrue(phenotyper.findGeneReport(g).orElseThrow(unfoundGene).isCalled(), g + " is not called"));
  }

  /**
   * Check to see if none of the given genes have been called by the matcher
   */
  private void assertNotCalledByMatcher(Phenotyper phenotyper, String... genes) {
    Arrays.stream(genes)
        .forEach(g -> assertFalse(phenotyper.findGeneReport(g).orElseThrow(unfoundGene).isCalled(), g + " is called"));
  }

  /**
   * Test the "print" calls for a gene that will display in the final report or in the phenotyper. This will check that
   * the call count matches and then check each individual call is present (can be 1 or more).
   * @param gene the gene to get diplotypes for
   * @param calls the expected display of the calls, 1 or more
   */
  private void assertDiplotypeDisplay(Phenotyper phenotyper, String gene, String... calls) {
    GeneReport geneReport = phenotyper.findGeneReport(gene).orElseThrow(unfoundGene);
    Collection<String> dips = geneReport.printDisplayCalls();
    assertEquals(calls.length, dips.size(),
        "Expected " + gene + " call count (" + calls.length + ") doesn't match actual call count (" +
            dips.size() + "): " + String.join(", ", dips));
    Arrays.stream(calls).forEach(c -> assertTrue(dips.contains(c), c + " not in " + gene + ":" + dips));
  }

  /**
   * Test the diplotype that will be used for looking up the recommendation. This will mostly match what's printed in
   * displays but will differ for particular genes
   * @param gene the gene to get diplotypes for
   * @param haplotypes the expected haplotypes names used for calling, specifying one will assume homozygous, otherwise specify two haplotype names
   */
  private void assertLookup(Phenotyper phenotyper, String gene, String... haplotypes) {
    Map<String,Integer> lookup = new HashMap<>();
    if (haplotypes.length == 1) {
      lookup.put(haplotypes[0], 2);
    } else if (haplotypes.length == 2) {
      if (haplotypes[0].equals(haplotypes[1])) {
        lookup.put(haplotypes[0], 2);
      } else {
        lookup.put(haplotypes[0], 1);
        lookup.put(haplotypes[1], 1);
      }
    } else {
      fail("Can only test on 1 or 2 haplotypes");
    }

    GeneReport geneReport = phenotyper.findGeneReport(gene).orElseThrow(unfoundGene);
    assertTrue(geneReport.isReportable());
    assertTrue(geneReport.getReporterDiplotypes().stream()
        .anyMatch(d -> d.makeLookupMap().equals(lookup)));
  }
}

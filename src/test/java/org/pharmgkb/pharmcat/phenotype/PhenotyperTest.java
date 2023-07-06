package org.pharmgkb.pharmcat.phenotype;

import java.net.URL;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.function.Supplier;
import com.google.common.collect.ImmutableList;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import org.pharmgkb.pharmcat.Env;
import org.pharmgkb.pharmcat.haplotype.ResultSerializer;
import org.pharmgkb.pharmcat.haplotype.model.GeneCall;
import org.pharmgkb.pharmcat.reporter.handlebars.ReportHelpers;
import org.pharmgkb.pharmcat.reporter.model.DataSource;
import org.pharmgkb.pharmcat.reporter.model.result.CallSource;
import org.pharmgkb.pharmcat.reporter.model.result.DiplotypeTest;
import org.pharmgkb.pharmcat.reporter.model.result.GeneReport;

import static org.junit.jupiter.api.Assertions.*;


/**
 * Unit test class for the {@link Phenotyper} class. Will run through synthetic data and assert genes are found and
 * diplotypes are called with certain results.
 */
class PhenotyperTest {
  private static final Supplier<RuntimeException> unfoundGene = () -> new RuntimeException("Gene report not found");
  private static Env s_env;

  @BeforeAll
  static void prepare() throws Exception {
    s_env = new Env();
  }

  @Test
  void testCyp2C19Het() throws Exception {
    Map<String, Collection<String>> warnings = new HashMap<>();
    warnings.put("chr10:94775453", ImmutableList.of("Test warning message"));
    warnings.put("chr10:94852914", ImmutableList.of("Test other message"));

    Phenotyper phenotyper = new Phenotyper(s_env,
        readMatchData("Cyp2C19Het.match.json"),
        OutsideCallParser.parse("CYP2D6\t*1/*3"), warnings);

    assertCalledByMatcher(phenotyper, "CYP2C19");
    assertReportable(phenotyper, "CYP2D6");
    assertNotCalledByMatcher(phenotyper, "CYP2C9");

    assertDiplotypeDisplay(phenotyper, "CYP2C19", "*1/*4");
    assertLookup(phenotyper, "CYP2C19", "*1", "*4");

    assertDiplotypeDisplay(phenotyper, "CYP2D6", "*1/*3");
    assertLookup(phenotyper, "CYP2D6", "*1", "*3");

    assertWarning(phenotyper, "rs72552267", "Test warning message");
    assertWarning(phenotyper, "rs55640102", "Test other message");
  }

  @Test
  void testCyp2D6Only() throws Exception {
    Phenotyper phenotyper = new Phenotyper(s_env,
        new ArrayList<>(),
        OutsideCallParser.parse("CYP2D6\t*1/*3"), null);

    assertReportable(phenotyper, "CYP2D6");

    // CYP2D6 is present and has an activity score
    assertTrue(phenotyper.findGeneReport(DataSource.CPIC, "CYP2D6").isPresent());
    phenotyper.findGeneReport(DataSource.CPIC, "CYP2D6").ifPresentOrElse(
        (geneReport) -> {
          assertTrue(geneReport.isCalled(), "CYP2D6 report should be present and called");
          assertEquals(CallSource.OUTSIDE, geneReport.getCallSource());
          assertEquals("1.0", geneReport.getRecommendationDiplotypes().first().getActivityScore());
        },
        () -> fail("CYP2D6 report should be present")
    );
    assertDiplotypeDisplay(phenotyper, "CYP2D6", "*1/*3");

    // CYP2C9 has a report but is not called and has no source
    assertTrue(phenotyper.findGeneReport(DataSource.CPIC, "CYP2C9").isPresent());
    phenotyper.findGeneReport(DataSource.CPIC, "CYP2C9").ifPresentOrElse(
        (geneReport) -> {
          assertEquals(CallSource.NONE, geneReport.getCallSource());
          assertFalse(geneReport.isCalled(), "CYP2C9 report should be present but not called");
        },
        () -> fail("CYP2C9 report should be present")
    );
  }

  @Test
  void testCyp2C19Hom() throws Exception {
    Phenotyper phenotyper = new Phenotyper(s_env,
        readMatchData("Cyp2C19s2s2.match.json"),
        new ArrayList<>(), null);

    assertCalledByMatcher(phenotyper, "CYP2C19");

    assertDiplotypeDisplay(phenotyper, "CYP2C19", "*2/*2");
  }

  @Test
  void testUGT1A1Phased() throws Exception {
    Phenotyper phenotyper = new Phenotyper(s_env,
        readMatchData("UGT1A1s1s60s80phased.match.json"),
        new ArrayList<>(), null);

    assertCalledByMatcher(phenotyper, "UGT1A1");

    assertDiplotypeDisplay(phenotyper, "UGT1A1", "*1/*80");
    assertLookup(phenotyper, "UGT1A1", "*1", "*80");
    assertTrue(phenotyper.findGeneReport(DataSource.CPIC, "UGT1A1").orElseThrow(unfoundGene).isPhased());
  }

  @Test
  void testUGT1A1Unphased() throws Exception {
    Phenotyper phenotyper = new Phenotyper(s_env,
        readMatchData("UGT1A1s1s60s80unphased.match.json"),
        new ArrayList<>(), null);

    assertCalledByMatcher(phenotyper, "UGT1A1");

    assertDiplotypeDisplay(phenotyper, "UGT1A1", "*1/*80");
    assertLookup(phenotyper, "UGT1A1", "*1", "*80");
    assertFalse(phenotyper.findGeneReport(DataSource.CPIC, "UGT1A1").orElseThrow(unfoundGene).isPhased());
  }

  @Test
  void testNUDT15() throws Exception {
    Phenotyper phenotyper = new Phenotyper(s_env,
        readMatchData("NUDT15ref.match.json"),
        new ArrayList<>(), null);

    assertCalledByMatcher(phenotyper, "NUDT15");

    assertDiplotypeDisplay(phenotyper, "NUDT15", "*1/*1");
    assertLookup(phenotyper, "NUDT15", "*1", "*1");
  }

  @Test
  void testNUDT15star3() throws Exception {
    Phenotyper phenotyper = new Phenotyper(s_env,
        readMatchData("NUDT15s3.match.json"),
        new ArrayList<>(), null);

    assertCalledByMatcher(phenotyper, "NUDT15");

    assertDiplotypeDisplay(phenotyper, "NUDT15", "*1/*3");
    assertLookup(phenotyper, "NUDT15", "*1", "*3");
  }


  //  Helper methods found below =======================================================================================

  private List<GeneCall> readMatchData(String testResourceFileName) throws Exception {
    URL testFileUrl = getClass().getResource(testResourceFileName);
    if (testFileUrl == null) {
      throw new RuntimeException("No test file found for " + testResourceFileName);
    }
    return new ResultSerializer().fromJson(Paths.get(testFileUrl.toURI())).getGeneCalls();
  }

  /**
   * Check to see if all the given genes have been called by the matcher
   */
  private void assertCalledByMatcher(Phenotyper phenotyper, String... genes) {
    Arrays.stream(genes)
        .forEach(g -> assertTrue(phenotyper.findGeneReport(DataSource.CPIC, g)
            .orElseThrow(unfoundGene).isCalled(), g + " is not called"));
  }

  private void assertReportable(Phenotyper phenotyper, String... genes) {
    Arrays.stream(genes)
        .forEach(g -> assertTrue(phenotyper.findGeneReport(DataSource.CPIC, g)
            .orElseThrow(unfoundGene).isReportable(), g + " is not reportable"));
  }

  /**
   * Check to see if none of the given genes have been called by the matcher
   */
  private void assertNotCalledByMatcher(Phenotyper phenotyper, String... genes) {
    Arrays.stream(genes)
        .forEach(g -> assertFalse(phenotyper.findGeneReport(DataSource.CPIC, g)
            .orElseThrow(unfoundGene).isCalled(), g + " is called"));
  }

  /**
   * Test the "print" calls for a gene that will display in the final report or in the phenotyper. This will check that
   * the call count matches and then check each individual call is present (can be 1 or more).
   * @param gene the gene to get diplotypes for
   * @param calls the expected display of the calls, 1 or more
   */
  private void assertDiplotypeDisplay(Phenotyper phenotyper, String gene, String... calls) {
    GeneReport geneReport = phenotyper.findGeneReport(DataSource.CPIC, gene).orElseThrow(unfoundGene);
    Collection<String> dips = ReportHelpers.amdGeneCalls(geneReport);
    assertEquals(calls.length, dips.size(),
        "Expected " + gene + " call count (" + calls.length + ") doesn't match actual call count (" +
            dips.size() + "): " + String.join(", ", dips));
    Arrays.stream(calls).forEach(c -> assertTrue(dips.contains(c), c + " not in " + gene + ":" + dips));
  }

  private void assertWarning(Phenotyper phenotyper, String rsid, String warning) {
    GeneReport geneReport = phenotyper.findGeneReport(DataSource.CPIC, "CYP2C19").orElseThrow(unfoundGene);
    assertTrue(geneReport.getVariantReports().stream()
        .anyMatch(v -> v.getDbSnpId() != null && v.getDbSnpId().equals(rsid) && v.getWarnings().contains(warning)));
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

    GeneReport geneReport = phenotyper.findGeneReport(DataSource.CPIC, gene).orElseThrow(unfoundGene);
    assertTrue(geneReport.isReportable());
    assertTrue(geneReport.getRecommendationDiplotypes().stream()
        .anyMatch(d -> DiplotypeTest.computeLookupMap(d).equals(lookup)));
  }
}

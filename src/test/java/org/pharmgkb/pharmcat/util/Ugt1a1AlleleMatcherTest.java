package org.pharmgkb.pharmcat.util;

import java.io.IOException;
import java.lang.invoke.MethodHandles;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Set;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import org.pharmgkb.common.util.PathUtils;
import org.pharmgkb.pharmcat.PharmCAT;
import org.pharmgkb.pharmcat.reporter.ReportContext;
import org.pharmgkb.pharmcat.reporter.model.result.GeneReport;

import static org.junit.jupiter.api.Assertions.*;


/**
 * Testing the specific functionality of the {@link Ugt1a1AlleleMatcher}
 *
 * @author Ryan Whaley
 */
class Ugt1a1AlleleMatcherTest {

  private static PharmCAT s_pharmcat;

  @BeforeAll
  static void prepare() throws IOException {
    s_pharmcat = new PharmCAT(Files.createTempDirectory(MethodHandles.lookup().lookupClass().getName()), null);
  }

  @Test
  void test() throws Exception {
    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/UGT1A1/s1s28s60s80unphased.vcf");
    s_pharmcat.execute(vcfFile, null, null);

    ReportContext context = s_pharmcat.getReporter().getContext();
    GeneReport geneReport = context.getGeneReport("UGT1A1");

    assertNotNull(geneReport);
    assertTrue(geneReport.isCalled());
    assertFalse(geneReport.isPhased());

    assertTrue(Ugt1a1AlleleMatcher.shouldBeUsedOn(geneReport));

    Set<String> lookupCalls = Ugt1a1AlleleMatcher.makeLookupCalls(geneReport);
    assertNotNull(lookupCalls);
    assertEquals(1, lookupCalls.size());
    assertTrue(lookupCalls.contains("*1/*80"));
  }

  @Test
  void testPhasedBalanced() throws Exception {
    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/UGT1A1/s6s60s80s28missingphased.vcf");
    s_pharmcat.execute(vcfFile, null, null);

    ReportContext context = s_pharmcat.getReporter().getContext();
    GeneReport geneReport = context.getGeneReport("UGT1A1");

    assertNotNull(geneReport);
    assertTrue(geneReport.isCalled());
    assertTrue(geneReport.isPhased());

    assertTrue(Ugt1a1AlleleMatcher.shouldBeUsedOn(geneReport));

    Set<String> lookupCalls = Ugt1a1AlleleMatcher.makeLookupCalls(geneReport);
    assertNotNull(lookupCalls);
    assertEquals(1, lookupCalls.size());
    assertTrue(lookupCalls.contains("*80/*80"));
  }

  @Test
  void testPhasedLopsided() throws Exception {
    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/UGT1A1/HG00436.vcf");
    s_pharmcat.execute(vcfFile, null, null);

    ReportContext context = s_pharmcat.getReporter().getContext();
    GeneReport geneReport = context.getGeneReport("UGT1A1");

    assertNotNull(geneReport);
    assertFalse(geneReport.isCalled());
    assertTrue(geneReport.isPhased());

    assertFalse(Ugt1a1AlleleMatcher.shouldBeUsedOn(geneReport));
  }
}

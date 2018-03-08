import java.io.FileWriter;
import java.io.IOException;
import java.lang.invoke.MethodHandles;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.function.Predicate;
import java.util.stream.Stream;
import org.junit.BeforeClass;
import org.junit.Test;
import org.pharmgkb.pharmcat.PharmCAT;
import org.pharmgkb.pharmcat.VcfTestUtils;
import org.pharmgkb.pharmcat.reporter.ReportContext;
import org.pharmgkb.pharmcat.reporter.model.VariantReport;
import org.pharmgkb.pharmcat.reporter.model.result.GeneReport;
import org.pharmgkb.pharmcat.reporter.model.result.GuidelineReport;

import static org.junit.Assert.*;


/**
 * Test the data generated from a full run of the PharmCAT matcher and reporter.
 *
 * @author Ryan Whaley
 */
public class PharmCATTest {

  private static final String sf_astrolabeOutput = "##Test Astrolabe output\n" +
      "#ROI_label\tdiplotype labels\tdiplotype activity\tdiplotype calling notes\tjaccard\tpart\tpValue\tROI notes\tspecial case\tnomenclature version\n" +
      "CYP2D6\tCYP2D6*1/CYP2D6*4\t?/?\t\t0.6\t0.75\tp: 0.0\t\t\tv1.9-2017_02_09\n";
  private static final String sf_diplotypesTemplate = "\nmatcher: %s\nreporter: %s\nprint (displayCalls): %s\nlookup: %s";
  private static PharmCAT s_pharmcat;
  private static Path s_tempAstroPath;
  private static ReportContext s_context;

  @BeforeClass
  public static void prepare() throws IOException {

    s_tempAstroPath = Files.createTempFile("astrolabe", ".tsv");
    try (FileWriter fw = new FileWriter(s_tempAstroPath.toFile())) {
      fw.write(sf_astrolabeOutput);
    }

    Path tempDirPath = Files.createTempDirectory(MethodHandles.lookup().lookupClass().getName());
    s_pharmcat = new PharmCAT(tempDirPath, null, null);
  }

  @Test
  public void testCyp2c19_1() throws Exception {
    generalTest("test.cyp2c19.s4s17het", new String[]{
            "cyp2c19/s4s17het.vcf"
        },
        false);


    testCalledGenes("CYP2C19");
    testCalls(DipType.PRINT,  "CYP2C19", "*1/*4B");

    testMatchedGroups("citalopram", 1);
    testMatchedGroups("ivacaftor", 0);
  }

  @Test
  public void testCyp2c19_astrolabe() throws Exception {
    generalTest("test.cyp2c19.s1s4b", new String[]{
        "cyp2c19/s4s17het.vcf"
        },
        true);

    testCalledGenes("CYP2C19", "CYP2D6");

    testCalls(DipType.PRINT, "CYP2D6", "*1/*4");
    testCalls(DipType.PRINT, "CYP2C19", "*1/*4B");

    assertTrue(s_context.getGeneReport("CYP2D6").isAstrolabeCall());
  }

  @Test
  public void testCftrRegInc() throws Exception {
    generalTest("test.cftr.reg_inc", new String[]{
            "CFTR/G542XF508del.vcf"
        },
        false);

    testCalledGenes("CFTR");
    testCalls(DipType.PRINT, "CFTR", "F508del(TCT)/G542X");
    testCalls(DipType.LOOKUP, "CFTR", "CFTR:F508del(TCT)/Other");

    assertTrue("Missing incidental allele", s_context.getGeneReports().stream().anyMatch(GeneReport::isIncidental));
  }

  @Test
  public void testCftrRefRef() throws Exception {
    generalTest("test.cftr.ref_ref", new String[]{
            "CFTR/refref.vcf"
        },
        false);

    testCalledGenes("CFTR");
    testCalls(DipType.PRINT, "CFTR", "No CPIC variants found");
    testCalls(DipType.LOOKUP, "CFTR", "CFTR:Other/Other");

    assertTrue("Should be no incidental alleles", s_context.getGeneReports().stream().noneMatch(GeneReport::isIncidental));
  }

  @Test
  public void testCftrF508() throws Exception {
    generalTest("test.cftr.refF508del", new String[]{
            "CFTR/refF508del.vcf"
        },
        false);

    testCalledGenes("CFTR");
    testCalls(DipType.PRINT, "CFTR", "F508del(TCT) (heterozygous)");
    testCalls(DipType.LOOKUP, "CFTR", "CFTR:F508del(TCT)/Other");

    assertTrue("Should be no incidental alleles", s_context.getGeneReports().stream().noneMatch(GeneReport::isIncidental));
  }

  @Test
  public void testCftrF508Hom() throws Exception {
    generalTest("test.cftr.F508delHom", new String[]{
            "CFTR/F508delF508del.vcf"
        },
        false);

    testCalledGenes("CFTR");
    testCalls(DipType.PRINT, "CFTR", "F508del(TCT)/F508del(TCT)");
    testCalls(DipType.LOOKUP, "CFTR", "CFTR:F508del(TCT)/F508del(TCT)");

    assertTrue("Should be no incidental alleles", s_context.getGeneReports().stream().noneMatch(GeneReport::isIncidental));
  }

  @Test
  public void testCftrF508HomCTT() throws Exception {
    generalTest("test.cftr.F508delHom_CTT", new String[]{
            "CFTR/F508delF508del_CTT.vcf"
        },
        false);

    testCalledGenes("CFTR");
    testCalls(DipType.PRINT, "CFTR", "F508del(CTT)/F508del(CTT)");
    testCalls(DipType.LOOKUP, "CFTR", "CFTR:F508del(CTT)/F508del(CTT)");

    assertTrue("Should be no incidental alleles", s_context.getGeneReports().stream().noneMatch(GeneReport::isIncidental));
  }

  /**
   * This is the same test case as {@link PharmCATTest#testCftrF508()} but the position lines have been sorted
   * lexigraphically. This should not affect the calling and should lead to same diplotype output.
   */
  @Test
  public void testCftrF508Sorted() throws Exception {
    generalTest("test.cftr.refF508del_sorted", new String[]{
            "CFTR/refF508del_sorted.vcf"
        },
        false);

    testCalledGenes("CFTR");
    testCalls(DipType.PRINT, "CFTR", "F508del(TCT) (heterozygous)");
    testCalls(DipType.LOOKUP, "CFTR", "CFTR:F508del(TCT)/Other");

    assertTrue("Should be no incidental alleles", s_context.getGeneReports().stream().noneMatch(GeneReport::isIncidental));
  }

  @Test
  public void testSlco1b1Test1() throws Exception {
    generalTest("test.slco1b1.17.21", new String[]{
            "SLCO1B1/s17s21.vcf"
        },
        false);

    testCalledGenes("SLCO1B1");
    testCalls(DipType.PRINT, "SLCO1B1", "*17/*21");
    testCalls(DipType.LOOKUP, "SLCO1B1", "SLCO1B1:*1A/*5");

    assertTrue("Should be no incidental alleles", s_context.getGeneReports().stream().noneMatch(GeneReport::isIncidental));
  }

  @Test
  public void testSlco1b1HomWild() throws Exception {
    generalTest("test.slco1b1.hom.wild", new String[]{
            "SLCO1B1/s1as1a.vcf"
        },
        false);

    testCalledGenes("SLCO1B1");
    testCalls(DipType.PRINT, "SLCO1B1", "*1A/*1A");
    testCalls(DipType.LOOKUP, "SLCO1B1", "SLCO1B1:*1A/*1A");

    assertTrue("Should be no incidental alleles", s_context.getGeneReports().stream().noneMatch(GeneReport::isIncidental));
  }

  @Test
  public void testSlco1b1HomVar() throws Exception {
    generalTest("test.slco1b1.hom.var", new String[]{
            "SLCO1B1/s5s15.vcf"
        },
        false);

    testCalledGenes("SLCO1B1");
    testCalls(DipType.PRINT, "SLCO1B1", "*5/*15");
    testCalls(DipType.LOOKUP, "SLCO1B1", "SLCO1B1:*5/*5");

    assertTrue("Should be no incidental alleles", s_context.getGeneReports().stream().noneMatch(GeneReport::isIncidental));
  }

  @Test
  public void testSlco1b1Test3() throws Exception {
    generalTest("test.slco1b1.1a.15", new String[]{
            "SLCO1B1/s1as15.vcf"
        },
        false);

    testCalledGenes("SLCO1B1");
    testCalls(DipType.PRINT, "SLCO1B1", "*1A/*15");
    testCalls(DipType.LOOKUP, "SLCO1B1", "SLCO1B1:*1A/*5");

    assertTrue("Should be no incidental alleles", s_context.getGeneReports().stream().noneMatch(GeneReport::isIncidental));
  }

  @Test
  public void testSlco1b1TestMissing() throws Exception {
    generalTest("test.slco1b1.missing", new String[]{
            "DPYD/s1s1.vcf",
            "UGT1A1/s1s1.vcf",
            "TPMT/s1s1.vcf"
        },
        false);

    testCalledGenes("DPYD", "UGT1A1", "TPMT");
    assertFalse(s_context.getGeneReport("SLCO1B1").isCalled());

    assertTrue("Should be no incidental alleles", s_context.getGeneReports().stream().noneMatch(GeneReport::isIncidental));
  }

  @Test
  public void testSlco1b1TestMulti() throws Exception {
    generalTest("test.slco1b1.multi", new String[]{
            "SLCO1B1/multi.vcf"
        },
        false);

    GeneReport geneReport = s_context.getGeneReport("SLCO1B1");
    assertNotNull(geneReport);
    assertFalse(geneReport.isCalled());

    testCalls(DipType.PRINT, "SLCO1B1", "rs4149056T/rs4149056C");
    testCalls(DipType.LOOKUP, "SLCO1B1", "SLCO1B1:*1A/*5");

    assertTrue("Should be no incidental alleles", s_context.getGeneReports().stream().noneMatch(GeneReport::isIncidental));
  }

  @Test
  public void testUgt1a1PhasedMulti() throws Exception {
    generalTest("test.ugt1a1.phased.multi", new String[]{
            "UGT1A1/s1s60s80phased.vcf"
        },
        false);

    testCalledGenes("UGT1A1");
    testCalls(DipType.PRINT, "UGT1A1", "*1/*60");
    testCalls(DipType.LOOKUP, "UGT1A1", "UGT1A1:*1/*60");

    assertTrue(s_context.getGeneReport("UGT1A1").isPhased());

    assertTrue("Should be no incidental alleles", s_context.getGeneReports().stream().noneMatch(GeneReport::isIncidental));
  }

  @Test
  public void testUgt1a1UnphasedMulti() throws Exception {
    generalTest("test.ugt1a1.unphased.multi", new String[]{
            "UGT1A1/s1s60s80unphased.vcf"
        },
        false);

    testCalledGenes("UGT1A1");
    testCalls(DipType.PRINT, "UGT1A1", "*1/*60");
    testCalls(DipType.LOOKUP, "UGT1A1", "UGT1A1:*1/*80");

    assertFalse(s_context.getGeneReport("UGT1A1").isPhased());

    assertTrue("Should be no incidental alleles", s_context.getGeneReports().stream().noneMatch(GeneReport::isIncidental));
  }

  @Test
  public void testTpmtStar1s() throws Exception {
    generalTest("test.tpmt.star1s", new String[]{
            "TPMT/s1ss1ss3.vcf"
        },
        false);

    testCalledGenes("TPMT");
    testCalls(DipType.PRINT, "TPMT", "*1/*3A");
    testCalls(DipType.LOOKUP, "TPMT", "TPMT:*1/*3A");

    GeneReport tpmtReport = s_context.getGeneReport("TPMT");
    assertEquals(30, tpmtReport.getVariantReports().size());
    assertEquals(1, tpmtReport.getVariantOfInterestReports().size());

    Predicate<VariantReport> singlePosition = r -> r.getDbSnpId() != null && r.getDbSnpId().equals("rs2842934");
    assertTrue(tpmtReport.getVariantOfInterestReports().stream().anyMatch(singlePosition));
    assertTrue(tpmtReport.getVariantOfInterestReports().stream().filter(singlePosition).allMatch(r -> r.getCall().equals("G|G")));

    assertEquals(0, tpmtReport.getHighlightedVariants().size());

    assertTrue("Should be no incidental alleles", s_context.getGeneReports().stream().noneMatch(GeneReport::isIncidental));
  }

  @Test
  public void testCombined() throws Exception {
    generalTest("test.combined", new String[]{
            "DPYD/s1s1.vcf",
            "UGT1A1/s1s1.vcf",
            "TPMT/s1s1.vcf",
            "CYP3A5/s1s7.vcf",
            "CFTR/refref.vcf",
            "CYP2C19/s2s2.vcf",
            "CYP2C9/s2s3.vcf",
            "SLCO1B1/s1as1a.vcf",
            "VKORC1/-1639A-1639A.vcf",
            "cyp4f2/s1s1.vcf",
            "IFNL3/rs12979860CC.vcf"
        },
        true);

    testCalledGenes("DPYD", "UGT1A1", "TPMT", "CYP3A5", "CFTR", "CYP2C19",
        "CYP2C9", "SLCO1B1", "VKORC1", "CYP4F2", "IFNL3", "CYP2D6");
    testCalls(DipType.PRINT, "TPMT", "*1/*1");
    testCalls(DipType.PRINT, "DPYD", "No CPIC decreased or no function variant with strong or moderate evidence found");
    testCalls(DipType.PRINT, "CYP2C19", "*2/*2");
    testCalls(DipType.LOOKUP, "TPMT", "TPMT:*1/*1");
    testCalls(DipType.PRINT, "CYP2D6", "*1/*4");

    assertTrue("Should be no incidental alleles", s_context.getGeneReports().stream().noneMatch(GeneReport::isIncidental));
  }


  /**
   * Runs the PharmCAT tool for the given example gene call data
   */
  private void generalTest(String name, String[] geneCalls, boolean includeAstrolabe) throws Exception {
    Path tempVcfPath = Files.createTempFile(name, ".vcf");
    try (FileWriter fw = new FileWriter(tempVcfPath.toFile())) {
      fw.write(VcfTestUtils.writeVcf(geneCalls));
    } catch (Exception ex) {
      ex.printStackTrace();
      throw ex;
    }

    Path astrolabePath = includeAstrolabe ? s_tempAstroPath : null;
    s_pharmcat.execute(tempVcfPath, astrolabePath, null);
    s_context = s_pharmcat.getReporter().getContext();

    assertEquals(14, s_context.getGeneReports().size());
    assertEquals(32, s_context.getGuidelineReports().size());
  }

  /**
   * Test the different types of diplotype calls that come out of the reporter
   * @param type what type of diplotype to test
   * @param gene the gene to get diplotypes for
   * @param calls the calls to match against
   */
  private void testCalls(DipType type, String gene, String... calls) {
    GeneReport geneReport = s_context.getGeneReport(gene);

    Collection<String> dips = type == DipType.PRINT ?
        geneReport.printDisplayCalls()
        : new ArrayList<>(geneReport.getDiplotypeLookupKeys());

    assertEquals(gene + " call count doesn't match", calls.length, dips.size());

    Arrays.stream(calls)
        .forEach(c -> assertTrue(c + " not in "+type+" for " + gene + ":" + dips + printDiagnostic(geneReport), dips.contains(c)));
  }

  private static String printDiagnostic(GeneReport geneReport) {
    return String.format(
        sf_diplotypesTemplate,
        geneReport.getMatcherDiplotypes().toString(),
        geneReport.getReporterDiplotypes().toString(),
        geneReport.printDisplayCalls(),
        geneReport.getDiplotypeLookupKeys()
    );
  }

  /**
   * Check to see if all the given genes have been called
   */
  private void testCalledGenes(String... genes) {
    assertTrue(genes != null && genes.length > 0);

    Arrays.stream(genes)
        .forEach(g -> assertTrue(g + " is not called", s_context.getGeneReport(g).isCalled()));
  }

  private void testMatchedGroups(String drugName, int count) {
    Stream<GuidelineReport> guidelineStream = s_context.getGuidelineReports().stream()
        .filter(r -> r.getRelatedDrugs().contains(drugName));

    if (count > 0) {
      assertTrue(
          drugName + " does not have matching group count of " + count,
          guidelineStream.allMatch(r -> r.getMatchingGroups().size() == count));
    }
    else {
      assertTrue(
          guidelineStream.allMatch(g -> g.getMatchingGroups() == null || g.getMatchingGroups().size() == 0));
    }
  }

  private enum DipType {
    PRINT, // the diplotype that is displayed to the end-user
    LOOKUP // the diplotype used to lookup annotations
  }
}

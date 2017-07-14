import java.io.FileWriter;
import java.io.IOException;
import java.lang.invoke.MethodHandles;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;
import org.junit.BeforeClass;
import org.junit.Test;
import org.pharmgkb.pharmcat.PharmCAT;
import org.pharmgkb.pharmcat.VcfTestUtils;
import org.pharmgkb.pharmcat.reporter.ReportContext;
import org.pharmgkb.pharmcat.reporter.model.result.Diplotype;
import org.pharmgkb.pharmcat.reporter.model.result.GeneReport;

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
  public void testSlco1b1Test1() throws Exception {
    generalTest("test.slco1b1.17.21", new String[]{
            "SLCO1B1/s17s21.vcf"
        },
        false);

    testCalledGenes("SLCO1B1");
    testCalls(DipType.PRINT, "SLCO1B1", "rs4149056T/rs4149056C");
    testCalls(DipType.LOOKUP, "SLCO1B1", "SLCO1B1:*1A/*5");

    assertTrue("Should be no incidental alleles", s_context.getGeneReports().stream().noneMatch(GeneReport::isIncidental));
  }

  @Test
  public void testSlco1b1Test2() throws Exception {
    generalTest("test.slco1b1.5.15", new String[]{
            "SLCO1B1/s5s15.vcf"
        },
        false);

    testCalledGenes("SLCO1B1");
    testCalls(DipType.PRINT, "SLCO1B1", "rs4149056C/rs4149056C");
    testCalls(DipType.LOOKUP, "SLCO1B1", "SLCO1B1:*5/*5");

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

    testCalledGenes("SLCO1B1");
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
    testCalls(DipType.PRINT, "UGT1A1", "*1/*60", "*1/*80");
    testCalls(DipType.LOOKUP, "UGT1A1", "UGT1A1:*1/*80");

    assertTrue(s_context.getGeneReport("UGT1A1").isPhased());

    assertTrue("Should be no incidental alleles", s_context.getGeneReports().stream().noneMatch(GeneReport::isIncidental));
  }

  @Test
  public void testUgt1a1UnphasedMulti() throws Exception {
    generalTest("test.ugt1a1.unphased.multi", new String[]{
            "UGT1A1/s1s60s80.vcf"
        },
        false);

    testCalledGenes("UGT1A1");
    testCalls(DipType.PRINT, "UGT1A1", "*1/*60", "*1/*80");
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

    assertTrue("Should be no incidental alleles", s_context.getGeneReports().stream().noneMatch(GeneReport::isIncidental));
  }


  /**
   * Runs the PharmCAT tool for the given example gene call data
   */
  private void generalTest(String name, String[] geneCalls, boolean includeAstrolabe) throws Exception {
    Path tempVcfPath = Files.createTempFile(name, ".vcf");
    try (FileWriter fw = new FileWriter(tempVcfPath.toFile())) {
      fw.write(VcfTestUtils.writeVcf(geneCalls));
    }

    Path astrolabePath = includeAstrolabe ? s_tempAstroPath : null;
    s_pharmcat.execute(tempVcfPath, astrolabePath, null);
    s_context = s_pharmcat.getReporter().getContext();

    assertEquals(14, s_context.getGeneReports().size());
    assertEquals(33, s_context.getGuidelineReports().size());
  }

  /**
   * Test the different types of diplotype calls that come out of the reporter
   * @param type what type of diplotype to test
   * @param gene the gene to get diplotypes for
   * @param calls the calls to match against
   */
  private void testCalls(DipType type, String gene, String... calls) {
    List<String> dips = type == DipType.PRINT ?
        s_context.getGeneReport(gene).getMatcherDiplotypes().stream().map(Diplotype::printDisplay).collect(Collectors.toList())
        : new ArrayList<>(s_context.getGeneReport(gene).getDiplotypeLookupKeys());

    assertEquals(gene + " call count doesn't match", calls.length, dips.size());

    Arrays.stream(calls)
        .forEach(c -> assertTrue(c + " not in "+type+" for " + gene + ":" + dips, dips.contains(c)));
  }

  /**
   * Check to see if all the given genes have been called
   */
  private void testCalledGenes(String... genes) {
    assertTrue(genes != null && genes.length > 0);

    Arrays.stream(genes)
        .forEach(g -> assertTrue(g + " is not called", s_context.getGeneReport(g).isCalled()));
  }

  private enum DipType {
    PRINT, // the diplotype that is displayed to the end-user
    LOOKUP // the diplotype used to lookup annotations
  }
}

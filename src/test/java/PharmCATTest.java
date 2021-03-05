import java.io.FileWriter;
import java.io.IOException;
import java.lang.invoke.MethodHandles;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import org.pharmgkb.pharmcat.ParseException;
import org.pharmgkb.pharmcat.PharmCAT;
import org.pharmgkb.pharmcat.VcfTestUtils;
import org.pharmgkb.pharmcat.reporter.ReportContext;
import org.pharmgkb.pharmcat.reporter.model.VariantReport;
import org.pharmgkb.pharmcat.reporter.model.result.DrugReport;
import org.pharmgkb.pharmcat.reporter.model.result.GeneReport;

import static org.junit.jupiter.api.Assertions.*;


/**
 * Test the data generated from a full run of the PharmCAT matcher and reporter.
 *
 * @author Ryan Whaley
 */
class PharmCATTest {

  private static final String sf_astrolabeOutput = "##Test Astrolabe output\n" +
      "#ROI_label\tdiplotype labels\tdiplotype activity\tdiplotype calling notes\tjaccard\tpart\tpValue\tROI notes\tspecial case\tnomenclature version\n" +
      "CYP2D6\tCYP2D6*1/CYP2D6*4\t?/?\t\t0.6\t0.75\tp: 0.0\t\t\tv1.9-2017_02_09\n";
  private static final String sf_diplotypesTemplate = "\nmatcher: %s\nreporter: %s\nprint (displayCalls): %s\nlookup: %s";
  private static PharmCAT s_pharmcat;
  private static Path s_tempAstroPath;
  private static ReportContext s_context;

  @BeforeAll
  static void prepare() throws IOException {

    s_tempAstroPath = Files.createTempFile("astrolabe", ".tsv");
    try (FileWriter fw = new FileWriter(s_tempAstroPath.toFile())) {
      fw.write(sf_astrolabeOutput);
    }

    Path tempDirPath = Files.createTempDirectory(MethodHandles.lookup().lookupClass().getName());
    s_pharmcat = new PharmCAT(tempDirPath, null);
  }

  @Test
  void testCyp2c19_1() throws Exception {
    generalTest("test.cyp2c19.s4s17het", new String[]{
            "cyp2c19/s4s17het.vcf"
        },
        null);


    testCalledGenes("CYP2C19");
    testCalls(DipType.PRINT,  "CYP2C19", "*1/*4");

    testMatchedGroups("citalopram", 1);
    testMatchedGroups("ivacaftor", 0);
  }

  @Test
  void testClomipramineCall() throws Exception {
    generalTest("test.clomipramine", new String[]{
            "cyp2c19/s2s2.vcf"
        },
        s_tempAstroPath);

    testCalledGenes("CYP2C19");
    testCalls(DipType.PRINT,  "CYP2C19", "*2/*2");

    testMatchedGroups("amitriptyline", 1);
    testMatchedGroups("clomipramine", 1);
    testMatchedGroups("desipramine", 1);
    testMatchedGroups("doxepin", 1);
    testMatchedGroups("imipramine", 1);
    testMatchedGroups("nortriptyline", 1);
    testMatchedGroups("trimipramine", 1);

    testMatchedGroups("clopidogrel", 1);

    testMatchedGroups("lansoprazole", 1);

    // voriconazole has 2 populations with recommendations so should have 2 matching groups
    testMatchedGroups("voriconazole", 2);
  }

  @Test
  void testCyp2c19noCall() throws Exception {
    generalTest("test.cyp2c19.noCall", new String[]{
            "cyp2c19/noCall.vcf"
        },
        null);


    assertFalse(s_context.getGeneReport("CYP2C19").isCalled());

    testMatchedGroups("citalopram", 0);
    testMatchedGroups("ivacaftor", 0);
  }

  @Test
  void testCyp2c19_astrolabe() throws Exception {
    generalTest("test.cyp2c19.s1s4b", new String[]{
        "cyp2c19/s4s17het.vcf"
        },
        s_tempAstroPath);

    testCalledGenes("CYP2C19", "CYP2D6");

    testCalls(DipType.PRINT, "CYP2D6", "*1/*4");
    testCalls(DipType.PRINT, "CYP2C19", "*1/*4");

    assertTrue(s_context.getGeneReport("CYP2D6").isOutsideCall());
  }

  @Test
  void testCftrRefRef() throws Exception {
    generalTest("test.cftr.ref_ref", new String[]{
            "cftr/refref.vcf"
        },
        null);

    testCalledGenes("CFTR");
    testCalls(DipType.PRINT, "CFTR", "No CPIC variants found");
    testCalls(DipType.LOOKUP, "CFTR", "CFTR:ivacaftor non-responsive CFTR sequence/ivacaftor non-responsive CFTR sequence");

    testMatchedGroups("ivacaftor", 1);
  }

  @Test
  void testCftrF508() throws Exception {
    generalTest("test.cftr.refF508del", new String[]{
            "cftr/refF508del.vcf"
        },
        null);

    testCalledGenes("CFTR");
    testCalls(DipType.PRINT, "CFTR", "No CPIC variants found");
    testCalls(DipType.LOOKUP, "CFTR", "CFTR:ivacaftor non-responsive CFTR sequence/ivacaftor non-responsive CFTR sequence");
  }

  @Test
  void testCftrF508HomCTT() throws Exception {
    generalTest("test.cftr.F508delHom_CTT", new String[]{
            "cftr/F508delF508del.vcf"
        },
        null);

    testCalledGenes("CFTR");
    testCalls(DipType.PRINT, "CFTR", "No CPIC variants found");
    testCalls(DipType.LOOKUP, "CFTR", "CFTR:ivacaftor non-responsive CFTR sequence/ivacaftor non-responsive CFTR sequence");
  }

  /**
   * This is the same test case as {@link PharmCATTest#testCftrF508()} but the position lines have been sorted
   * lexigraphically. This should not affect the calling and should lead to same diplotype output.
   */
  @Test
  void testCftrF508Sorted() throws Exception {
    generalTest("test.cftr.refF508del_sorted", new String[]{
            "cftr/refF508del_sorted.vcf"
        },
        null);

    testCalledGenes("CFTR");
    testCalls(DipType.PRINT, "CFTR", "No CPIC variants found");
    testCalls(DipType.LOOKUP, "CFTR", "CFTR:ivacaftor non-responsive CFTR sequence/ivacaftor non-responsive CFTR sequence");
  }

  @Test
  void testSlco1b1Test1() throws Exception {
    generalTest("test.slco1b1.17.21", new String[]{
            "SLCO1B1/s17s21.vcf"
        },
        null);

    testCalledGenes("SLCO1B1");
    testCalls(DipType.PRINT, "SLCO1B1", "*17/*21");
    testCalls(DipType.LOOKUP, "SLCO1B1", "SLCO1B1:*17/*21");
  }

  @Test
  void testSlco1b1HomWild() throws Exception {
    generalTest("test.slco1b1.hom.wild", new String[]{
            "SLCO1B1/s1as1a.vcf"
        },
        null);

    testCalledGenes("SLCO1B1");
    testCalls(DipType.PRINT, "SLCO1B1", "*1A/*1A");
    testCalls(DipType.LOOKUP, "SLCO1B1", "SLCO1B1:*1A/*1A");
  }

  @Test
  void testSlco1b1HomVar() throws Exception {
    generalTest("test.slco1b1.hom.var", new String[]{
            "SLCO1B1/s5s15.vcf"
        },
        null);

    testCalledGenes("SLCO1B1");
    testCalls(DipType.PRINT, "SLCO1B1", "*5/*15");
    testCalls(DipType.LOOKUP, "SLCO1B1", "SLCO1B1:*5/*15");
  }

  @Test
  void testSlco1b1Test3() throws Exception {
    generalTest("test.slco1b1.1a.15", new String[]{
            "SLCO1B1/s1as15.vcf"
        },
        null);

    testCalledGenes("SLCO1B1");
    testCalls(DipType.PRINT, "SLCO1B1", "*1A/*15");
    testCalls(DipType.LOOKUP, "SLCO1B1", "SLCO1B1:*1A/*15");
  }

  @Test
  void testSlco1b1TestMissing() throws Exception {
    generalTest("test.slco1b1.missing", new String[]{
            "DPYD/s1s1.vcf",
            "TPMT/s1s1.vcf"
        },
        null);

    testCalledGenes("DPYD", "TPMT");
    assertFalse(s_context.getGeneReport("SLCO1B1").isCalled());
  }

  @Test
  void testDpydS1S2B() throws Exception {
    generalTest("test.slco1b1.missing", new String[]{
            "DPYD/s1s2b.vcf"
        },
        null);

    testCalledGenes("DPYD");
    assertTrue(s_context.getGeneReport("DPYD").isCalled());
    testCalls(DipType.PRINT, "DPYD", "c.1627A>G (*5)/c.1905+1G>A (*2A)");
    testCalls(DipType.LOOKUP, "DPYD", "DPYD:c.1627A>G (*5)/c.1905+1G>A (*2A)");

    testMatchedGroups("fluorouracil", 1);
    testMatchedGroups("capecitabine", 1);
  }

  @Test
  void testSlco1b1TestMulti() throws Exception {
    generalTest("test.slco1b1.multi", new String[]{
            "SLCO1B1/multi.vcf"
        },
        null);

    GeneReport geneReport = s_context.getGeneReport("SLCO1B1");
    assertNotNull(geneReport);
    assertFalse(geneReport.isCalled());

    testCalls(DipType.PRINT, "SLCO1B1", "rs4149056T/rs4149056C");
    testCalls(DipType.LOOKUP, "SLCO1B1", "SLCO1B1:*1A/*5");

    testMatchedGroups("simvastatin", 1);
  }

  @Test
  void testUgt1a1PhasedMulti() throws Exception {
    generalTest("test.ugt1a1.phased.multi", new String[]{
            "UGT1A1/s1s60s80phased.vcf"
        },
        null);

    testCalledGenes("UGT1A1");
    testCalls(DipType.PRINT, "UGT1A1", "*1/*80");
    testCalls(DipType.LOOKUP, "UGT1A1", "UGT1A1:*1/*80");

    assertTrue(s_context.getGeneReport("UGT1A1").isPhased());
  }

  @Test
  void testUgt1a1UnphasedMulti() throws Exception {
    generalTest("test.ugt1a1.unphased.multi", new String[]{
            "UGT1A1/s1s60s80unphased.vcf"
        },
        null);

    testCalledGenes("UGT1A1");
    testCalls(DipType.PRINT, "UGT1A1", "*80 (heterozygous)");
    testCalls(DipType.LOOKUP, "UGT1A1", "UGT1A1:*1/*1");

    assertFalse(s_context.getGeneReport("UGT1A1").isPhased());
  }

  @Test
  void testUgt1a1S1S28S60S80() throws Exception {
    generalTest("test.ugt1a1.s1s28s60s80unphased", new String[]{
            "UGT1A1/s1s28s60s80unphased.vcf"
        },
        null);

    testCalledGenes("UGT1A1");
    testCalls(DipType.PRINT, "UGT1A1", "*80+*28 (heterozygous)", "*28 (heterozygous)", "*80 (heterozygous)");
    testCalls(DipType.LOOKUP, "UGT1A1", "UGT1A1:*1/*80");

    assertFalse(s_context.getGeneReport("UGT1A1").isPhased());
  }

  @Test
  void testUgt1a1S28S37() throws Exception {
    generalTest("test.ugt1a1.s28s37", new String[]{
            "UGT1A1/s28s37.vcf"
        },
        null);

    testCalledGenes("UGT1A1");
    testCalls(DipType.PRINT, "UGT1A1", "*28 (heterozygous)", "*37 (heterozygous)");
    testCalls(DipType.LOOKUP, "UGT1A1", "UGT1A1:*80/*80");

    assertFalse(s_context.getGeneReport("UGT1A1").isPhased());
  }

  @Test
  void testUgt1a1s28s80phased() throws Exception {
    generalTest("test.ugt1a1.s28s80phased", new String[]{
            "UGT1A1/s28s80phased.vcf"
        },
        null);

    testCalledGenes("UGT1A1");
    testCalls(DipType.PRINT, "UGT1A1", "*1/*28+*80");
    testCalls(DipType.LOOKUP, "UGT1A1", "UGT1A1:*1/*80");

    assertTrue(s_context.getGeneReport("UGT1A1").isPhased());
  }

  @Test
  void testUgt1a1s28s80s6s60phased() throws Exception {
    generalTest("test.ugt1a1.s28s80s6s60phased", new String[]{
            "UGT1A1/s28s80s6s60phased.vcf"
        },
        null);

    testCalledGenes("UGT1A1");
    testCalls(DipType.PRINT, "UGT1A1", "*6/*28+*80");
    testCalls(DipType.LOOKUP, "UGT1A1", "UGT1A1:*80/*80");

    assertTrue(s_context.getGeneReport("UGT1A1").isPhased());
  }

  @Test
  void testUgt1a1s28s80s6s60unphased() throws Exception {
    generalTest("test.ugt1a1.s28s80s6s60unphased", new String[]{
            "UGT1A1/s28s80s6s60unphased.vcf"
        },
        null);

    testCalledGenes("UGT1A1");
    testCalls(DipType.PRINT, "UGT1A1", "*80+*28 (heterozygous)","*28 (heterozygous)","*6 (heterozygous)","*80 (heterozygous)");
    testCalls(DipType.LOOKUP, "UGT1A1", "UGT1A1:*80/*80");

    assertFalse(s_context.getGeneReport("UGT1A1").isPhased());
  }

  @Test
  void testUgt1a1s28s80unphased() throws Exception {
    generalTest("test.ugt1a1.s28s80unphased", new String[]{
            "UGT1A1/s28s80unphased.vcf"
        },
        null);

    testCalledGenes("UGT1A1");
    testCalls(DipType.PRINT, "UGT1A1", "*80+*28 (heterozygous)", "*28 (heterozygous)", "*80 (heterozygous)");
    testCalls(DipType.LOOKUP, "UGT1A1", "UGT1A1:*1/*80");

    assertFalse(s_context.getGeneReport("UGT1A1").isPhased());
  }

  @Test
  void testUgt1a1s6s6() throws Exception {
    generalTest("test.ugt1a1.s6s6", new String[]{
            "UGT1A1/s6s6.vcf"
        },
        null);

    testCalledGenes("UGT1A1");
    testCalls(DipType.PRINT, "UGT1A1", "*6/*6");
    testCalls(DipType.LOOKUP, "UGT1A1", "UGT1A1:*6/*6");

    assertTrue(s_context.getGeneReport("UGT1A1").isPhased());
  }

  @Test
  void testUgt1a1s6s60s80s28MissingPhased() throws Exception {
    generalTest("test.ugt1a1.s6s60s80s28MissingPhased", new String[]{
            "UGT1A1/s6s60s80s28missingphased.vcf"
        },
        null);

    testCalledGenes("UGT1A1");
    testCalls(DipType.PRINT, "UGT1A1", "*6/*28+*37+*80");
    testCalls(DipType.LOOKUP, "UGT1A1", "UGT1A1:*80/*80");

    assertTrue(s_context.getGeneReport("UGT1A1").isPhased());
  }

  @Test
  void testUgt1a1s6s60s80s28MissingUnphased() throws Exception {
    generalTest("test.ugt1a1.s6s60s80s28MissingUnphased", new String[]{
            "UGT1A1/s6s60s80s28missingunphased.vcf"
        },
        null);

    testCalledGenes("UGT1A1");
    testCalls(DipType.PRINT, "UGT1A1", "*6 (heterozygous)","*80 (heterozygous)","*80+*28 (heterozygous)","*80+*37 (heterozygous)");
    testCalls(DipType.LOOKUP, "UGT1A1", "UGT1A1:*80/*80");

    assertFalse(s_context.getGeneReport("UGT1A1").isPhased());
  }

  @Test
  void testUgt1a1s80s28missing() throws Exception {
    generalTest("test.ugt1a1.s80s28missing", new String[]{
            "UGT1A1/s80s28missing.vcf"
        },
        null);

    testCalledGenes("UGT1A1");
    testCalls(DipType.PRINT, "UGT1A1", "*80 (heterozygous)", "*80+*28 (heterozygous)", "*80+*37 (heterozygous)");
    testCalls(DipType.LOOKUP, "UGT1A1", "UGT1A1:*1/*80");

    assertFalse(s_context.getGeneReport("UGT1A1").isPhased());
  }

  @Test
  void testUgt1a1na12717() throws Exception {
    generalTest("test.ugt1a1.na12717", new String[]{
            "UGT1A1/NA12717_UGT1A1.vcf"
        },
        null);

    testCalledGenes("UGT1A1");
    testCalls(DipType.PRINT, "UGT1A1", "*80+*28 (heterozygous)", "*80 (heterozygous)", "*80 (homozygous)", "*28 (heterozygous)");
    testCalls(DipType.LOOKUP, "UGT1A1", "UGT1A1:*1/*80");

    assertFalse(s_context.getGeneReport("UGT1A1").isPhased());
  }

  @Test
  void testUgt1a1na18868() throws Exception {
    generalTest("test.ugt1a1.na18868", new String[]{
            "UGT1A1/NA18868_UGT1A1.vcf"
        },
        null);

    testCalledGenes("UGT1A1");
    testCalls(DipType.PRINT, "UGT1A1", "*80+*28 (heterozygous)", "*80 (heterozygous)", "*80 (homozygous)", "*28 (heterozygous)");
    testCalls(DipType.LOOKUP, "UGT1A1", "UGT1A1:*1/*80");

    assertFalse(s_context.getGeneReport("UGT1A1").isPhased());
  }

  @Test
  void testUgt1a1na19785() throws Exception {
    generalTest("test.ugt1a1.na19785", new String[]{
            "UGT1A1/NA19785_UGT1A1.vcf"
        },
        null);

    testCalledGenes("UGT1A1");
    testCalls(DipType.PRINT, "UGT1A1", "*80+*28 (heterozygous)", "*80 (heterozygous)", "*80 (homozygous)", "*28 (heterozygous)");
    testCalls(DipType.LOOKUP, "UGT1A1", "UGT1A1:*1/*80");

    assertFalse(s_context.getGeneReport("UGT1A1").isPhased());
  }

  @Test
  void testUgt1a1s28homMissing() throws Exception {
    generalTest("test.ugt1a1.s28s28unphaseds60s80miss", new String[]{
            "UGT1A1/s28s28unphaseds60s80miss.vcf"
        },
        null);

    testCalledGenes("UGT1A1");
    testCalls(DipType.PRINT, "UGT1A1", "*28+*80/*28+*80");
    testCalls(DipType.LOOKUP, "UGT1A1", "UGT1A1:*80/*80");

    // sample is effectively phased since all positions homozygous
    assertTrue(s_context.getGeneReport("UGT1A1").isPhased());
  }

  /**
   * Example of a UGT1A1 report that looks odd. The matcher calls *28/*60 and *60/*60 which then gets translated to the
   * print diplotypes listed below. It looks odd to have both "*60 (homozygous)" and "*60 (heterozygous)". This happens
   * because diplotype calls get translated to the zygosity format one at a time and can't "look ahead" to other
   * matches to check for homozygous calls that could obviate a heterozygous call display.
   *
   * Leaving this here for now but could be addressed in a future release.
   */
  @Test
  void testUgt1a1s28s60Hom() throws Exception {
    generalTest("test.ugt1a1.s28s60hom", new String[]{
            "UGT1A1/s28s60hom.vcf"
        },
        null);

    testCalledGenes("UGT1A1");
    testCalls(DipType.PRINT, "UGT1A1", "*28 (heterozygous)");
    testCalls(DipType.LOOKUP, "UGT1A1", "UGT1A1:*1/*80");

    // sample is effectively phased since all positions homozygous
    assertFalse(s_context.getGeneReport("UGT1A1").isPhased());
  }

  @Test
  void testUgt1a1s27s28unphaseds80s60missing() throws Exception {
    generalTest("test.ugt1a1.s27s28unphaseds80s60missing", new String[]{
            "UGT1A1/s27s28unphaseds80s60missing.vcf"
        },
        null);

    testCalledGenes("UGT1A1");
    testCalls(DipType.PRINT, "UGT1A1", "*27 (heterozygous)", "*28 (heterozygous)", "*80+*28 (heterozygous)");
    testCalls(DipType.LOOKUP, "UGT1A1", "UGT1A1:*80/*80");

    assertFalse(s_context.getGeneReport("UGT1A1").isPhased());
  }

  @Test
  void testUgt1a1s28hets60homounphaseds80missing() throws Exception {
    generalTest("test.ugt1a1.s28hets60homounphaseds80missing", new String[]{
            "UGT1A1/s28hets60homounphaseds80missing.vcf"
        },
        null);

    testCalledGenes("UGT1A1");
    testCalls(DipType.PRINT, "UGT1A1", "*28 (heterozygous)", "*80+*28 (heterozygous)");
    testCalls(DipType.LOOKUP, "UGT1A1", "UGT1A1:*1/*80");

    assertFalse(s_context.getGeneReport("UGT1A1").isPhased());
  }

  @Test
  void testUgt1a1HG00436() throws Exception {
    generalTest("test.ugt1a1.HG00436", new String[]{
            "UGT1A1/HG00436.vcf"
        },
        null);

    testCalledGenes("UGT1A1");
    testCalls(DipType.PRINT, "UGT1A1", "*1/*27+*28+*80");
    testCalls(DipType.LOOKUP, "UGT1A1", "UGT1A1:*1/*80");

    GeneReport geneReport = s_context.getGeneReport("UGT1A1");
    assertTrue(geneReport.isPhased());
  }

  @Test
  void testUgt1a1s1s80s27s60s28missingphased() throws Exception {
    generalTest("test.ugt1a1.s1s80s27s60s28missingphased", new String[]{
            "UGT1A1/s1s80s27s60s28missingphased.vcf"
        },
        null);

    testCalledGenes("UGT1A1");
    testCalls(DipType.PRINT, "UGT1A1", "*1/*27+*28+*37+*80");
    testCalls(DipType.LOOKUP, "UGT1A1", "UGT1A1:*1/*80");

    GeneReport geneReport = s_context.getGeneReport("UGT1A1");
    assertTrue(geneReport.isPhased());
  }

  @Test
  void testUgt1a1s1s60s80s6phased() throws Exception {
    generalTest("test.ugt1a1.s1s60s80s6phased", new String[]{
            "UGT1A1/s1s60s80s6phased.vcf"
        },
        null);

    testCalledGenes("UGT1A1");
    testCalls(DipType.PRINT, "UGT1A1", "*1/*6+*80");
    testCalls(DipType.LOOKUP, "UGT1A1", "UGT1A1:*1/*80");

    GeneReport geneReport = s_context.getGeneReport("UGT1A1");
    assertTrue(geneReport.isPhased());
  }

  @Test
  void testUgt1a1s1s60s80s28s6phased() throws Exception {
    generalTest("test.ugt1a1.s1s60s80s28s6phased", new String[]{
            "UGT1A1/s1s60s80s28s6phased.vcf"
        },
        null);

    testCalledGenes("UGT1A1");
    testCalls(DipType.PRINT, "UGT1A1", "*1/*6+*28+*80");
    testCalls(DipType.LOOKUP, "UGT1A1", "UGT1A1:*1/*80");

    GeneReport geneReport = s_context.getGeneReport("UGT1A1");
    assertTrue(geneReport.isPhased());
  }

  @Test
  void testUgt1a1s1s37s80s60phased() throws Exception {
    generalTest("test.ugt1a1.s1s37s80s60phased", new String[]{
            "UGT1A1/s1s37s80s60phased.vcf"
        },
        null);

    testCalledGenes("UGT1A1");
    testCalls(DipType.PRINT, "UGT1A1", "*1/*37+*80");
    testCalls(DipType.LOOKUP, "UGT1A1", "UGT1A1:*1/*80");

    GeneReport geneReport = s_context.getGeneReport("UGT1A1");
    assertTrue(geneReport.isPhased());
  }

  @Test
  void testCyp3a5Missing3Message() throws Exception {
    String gene = "CYP3A5";

    generalTest("test.cyp3a5.s3missing", new String[]{
            "cyp3a5/s1s1rs776746missing.vcf"
        },
        null);

    testCalledGenes(gene);
    testCalls(DipType.PRINT, gene, "*1/*1");
    testCalls(DipType.LOOKUP, gene, "CYP3A5:*1/*1");

    // rs776746 should be missing from this report
    assertNotNull(s_context.getGeneReport(gene).getVariantReports());
    assertTrue(s_context.getGeneReport(gene).getVariantReports().stream().anyMatch(v -> v.isMissing() && v.getDbSnpId().equals("rs776746")));

    // the guideline should have a matching message
    assertTrue(s_context.getDrugReports().stream()
        .filter(r -> r.getRelatedDrugs().contains("tacrolimus"))
        .allMatch(r -> r.getMessages().size() > 0));

    assertTrue(s_context.getGeneReport(gene).isPhased());
  }

  @Test
  void testCyp3a5MissingRS776746() throws Exception {
    generalTest("test.cyp3a5.missingRs776746", new String[]{
            "cyp3a5/s1s1rs776746missing.vcf"
        },
        null);
    
    testCalledGenes("CYP3A5");
    testCalls(DipType.PRINT, "CYP3A5", "*1/*1");
  }

  @Test
  void testCyp3a5v1() throws Exception {
    generalTest("test.cyp3a5.s1s3rs776746rs55965422het", new String[]{
            "cyp3a5/s1s3rs776746rs55965422het.vcf"
        },
        null);
    
    testCalledGenes("CYP3A5");
    testCalls(DipType.PRINT, "CYP3A5", "*1/*3");
  }

  @Test
  void testCyp3a5v2() throws Exception {
    generalTest("test.cyp3a5.s1s3rs776746rs55965422rs28383479het", new String[]{
            "cyp3a5/s1s3rs776746rs55965422rs28383479het.vcf"
        },
        null);
    
    testCalledGenes("CYP3A5");
    testCalls(DipType.PRINT, "CYP3A5", "*1/*3");
  }

  @Test
  void testCyp3a5v3() throws Exception {
    generalTest("test.cyp3a5.s3s3rs55965422het", new String[]{
            "cyp3a5/s3s3rs55965422het.vcf"
        },
        null);
    
    testCalledGenes("CYP3A5");
    testCalls(DipType.PRINT, "CYP3A5", "*3/*3");
  }

  @Test
  void testCyp3a5v4() throws Exception {
    generalTest("test.cyp3a5.s3s5-homozygous", new String[]{
            "cyp3a5/s3s5-homozygous.vcf"
        },
        null);
    
    testCalledGenes("CYP3A5");
    testCalls(DipType.PRINT, "CYP3A5", "*3/*5");
  }

  @Test
  void testCyp3a5v5() throws Exception {
    generalTest("test.cyp3a5.s1s3rs776746rs28383479het", new String[]{
            "cyp3a5/s1s3rs776746rs28383479het.vcf"
        },
        null);
    
    testCalledGenes("CYP3A5");
    testCalls(DipType.PRINT, "CYP3A5", "*1/*3");
  }

  @Test
  void testTpmtStar1s() throws Exception {
    generalTest("test.tpmt.star1s", new String[]{
            "TPMT/s1ss1ss3.vcf"
        },
        null);

    testCalledGenes("TPMT");
    testCalls(DipType.PRINT, "TPMT", "*1/*3A");
    testCalls(DipType.LOOKUP, "TPMT", "TPMT:*1/*3A");

    GeneReport tpmtReport = s_context.getGeneReport("TPMT");
    assertEquals(43, tpmtReport.getVariantReports().size());

    assertEquals(0, tpmtReport.getHighlightedVariants().size());
  }

  @Test
  void testTpmtS15OffData() throws Exception {
    generalTest("test.tpmt.s15offdata", new String[] {
            "TPMT/s15offdata.vcf"
        },
        null);

    testNotCalledGenes("TPMT");
    GeneReport report = s_context.getGeneReport("TPMT");
    assertTrue(report.getVariantReports().stream().filter(r -> r.getPosition() == 18133890).allMatch(VariantReport::isMismatch));
  }


  @Test
  void testCyp2c9star61() throws Exception {
    generalTest("test.cyp2c9.s1s61", new String[] {
            "cyp2c9/s1s61.vcf"
        },
        null);

    testCalledGenes("CYP2C9");
    testCalls(DipType.PRINT, "CYP2C9", "*1/*61");
    testCalls(DipType.LOOKUP, "CYP2C9", "CYP2C9:*1/*61");

    testMatchedGroups("lornoxicam", 1);
  }

  @Test
  void testCyp2c9star1Hom() throws Exception {
    generalTest("test.cyp2c9.s1s1", new String[] {
            "cyp2c9/s1s1.vcf"
        },
        null);

    testCalledGenes("CYP2C9");
    testCalls(DipType.PRINT, "CYP2C9", "*1/*1");
    testCalls(DipType.LOOKUP, "CYP2C9", "CYP2C9:*1/*1");

    testMatchedGroups("celecoxib", 1);
    testMatchedGroups("ibuprofen", 1);
    testMatchedGroups("lornoxicam", 1);
  }


  @Test
  void testCombined() throws Exception {
    generalTest("test.combined", new String[]{
            "DPYD/s1s1.vcf",
            "UGT1A1/s1s1.vcf",
            "TPMT/s1s1.vcf",
            "cyp3a5/s1s7.vcf",
            "cftr/refref.vcf",
            "cyp2c19/s2s2.vcf",
            "cyp2c9/s2s3.vcf",
            "SLCO1B1/s1as1a.vcf",
            "VKORC1/-1639A-1639A.vcf",
            "cyp4f2/s1s1.vcf",
            "IFNL3/rs12979860CC.vcf"
        },
        s_tempAstroPath);

    testCalledGenes("DPYD", "UGT1A1", "TPMT", "CYP3A5", "CFTR", "CYP2C19",
        "CYP2C9", "SLCO1B1", "VKORC1", "CYP4F2", "IFNL3", "CYP2D6");
    testCalls(DipType.PRINT, "TPMT", "*1/*1");
    testCalls(DipType.PRINT, "DPYD", "Reference/Reference");
    testCalls(DipType.PRINT, "CYP2C19", "*2/*2");
    testCalls(DipType.LOOKUP, "TPMT", "TPMT:*1/*1");
    testCalls(DipType.PRINT, "CYP2D6", "*1/*4");
    testCalls(DipType.PRINT, "UGT1A1", "*1/*1");
  }

  @Test
  void testBadOutsideData() throws Exception {
    Path badOutsideDataPath = Files.createTempFile("astrolabe", ".tsv");
    try (FileWriter fw = new FileWriter(badOutsideDataPath.toFile())) {
      fw.write("CYP2D6\t*1/*2\nCYP2D6\t*3/*4");
    }

    try {
      generalTest("test.badOutsideData", new String[]{
          "cyp2c19/s2s2.vcf",
          "cyp2c9/s2s3.vcf",
      }, badOutsideDataPath);
      fail("Should have failed due to a duplicate gene definition in outside call");
    }
    catch (ParseException ex) {
      // we want this to fail so ignore handling the exception
    }
  }

  @Test
  void testCallerCollision() throws Exception {
    Path outsidePath = Files.createTempFile("astrolabe", ".tsv");
    try (FileWriter fw = new FileWriter(outsidePath.toFile())) {
      fw.write("CYP2C19\t*2/*2");
    }

    try {
      generalTest("test.badOutsideData", new String[]{
          "cyp2c19/s2s2.vcf",
      }, outsidePath);
      fail("Should have failed due to a duplicate gene definition between matcher and outside caller");
    }
    catch (ParseException ex) {
      // we want this to fail so ignore handling the exception
    }
  }


  /**
   * Runs the PharmCAT tool for the given example gene call data
   */
  private void generalTest(String name, String[] geneCalls, Path outsideCallPath) throws Exception {
    Path tempVcfPath = Files.createTempFile(name, ".vcf");
    try (FileWriter fw = new FileWriter(tempVcfPath.toFile())) {
      fw.write(VcfTestUtils.writeVcf(geneCalls));
    } catch (Exception ex) {
      ex.printStackTrace();
      throw ex;
    }

    s_pharmcat.execute(tempVcfPath, outsideCallPath, null);
    s_context = s_pharmcat.getReporter().getContext();

    assertEquals(16, s_context.getGeneReports().size());
    assertEquals(49, s_context.getDrugReports().size());
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

    assertEquals(calls.length, dips.size(), gene + " call count doesn't match " + String.join(";", dips));

    Arrays.stream(calls)
        .forEach(c -> assertTrue(dips.contains(c),
            c + " not in "+type+" for " + gene + ":" + dips + printDiagnostic(geneReport)));
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
        .forEach(g -> assertTrue(s_context.getGeneReport(g).isCalled(), g + " is not called"));
  }

  /**
   * Check to see if none ofthe given genes have been called
   */
  private void testNotCalledGenes(String... genes) {
    assertTrue(genes != null && genes.length > 0);

    Arrays.stream(genes)
        .forEach(g -> assertFalse(s_context.getGeneReport(g).isCalled(), g + " is called"));
  }

  /**
   * Check to see if there is a matching recommendation for the given drug name
   * @param drugName a drug name that has recommendations
   * @param expectedCount the number of matching recommendations you expect
   */
  private void testMatchedGroups(String drugName, int expectedCount) {
    DrugReport guideline = s_context.getDrugReports().stream()
        .filter(r -> r.getRelatedDrugs().contains(drugName))
        .findFirst().orElseThrow(() -> new RuntimeException("No guideline found for " + drugName));

    assertEquals(expectedCount, guideline.getMatchingRecommendations().size(),
        drugName + " does not have matching recommendation count of " + expectedCount);
  }

  private enum DipType {
    PRINT, // the diplotype that is displayed to the end-user
    LOOKUP // the diplotype used to lookup annotations
  }
}

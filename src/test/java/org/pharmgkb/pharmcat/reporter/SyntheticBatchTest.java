package org.pharmgkb.pharmcat.reporter;

import java.io.FileWriter;
import java.io.IOException;
import java.lang.invoke.MethodHandles;
import java.nio.file.Files;
import java.nio.file.Path;
import org.pharmgkb.common.util.CliHelper;
import org.pharmgkb.common.util.PathUtils;
import org.pharmgkb.pharmcat.PharmCAT;
import org.pharmgkb.pharmcat.VcfTestUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


/**
 * This test generates a collection of synthetic sample VCF files from our test resources and then runs the matcher and
 * reporter on those samples, writing the output to a desired directory.
 *
 * This is especially useful to see how different allele calls will be called and displayed in the final output report.
 *
 * @author Ryan Whaley
 */
class SyntheticBatchTest {
  private static final Logger sf_logger = LoggerFactory.getLogger(MethodHandles.lookup().lookupClass());
  private static final Path sf_outsideCYP2D6File
      = PathUtils.getPathToResource("org/pharmgkb/pharmcat/reporter/outside_CYP2D6.tsv");
  private static final Path sf_outsideCYP2D6G6PDFile
      = PathUtils.getPathToResource("org/pharmgkb/pharmcat/reporter/outside_CYP2D6_G6PD.tsv");

  private final PharmCAT f_pharmcat;
  private final Path f_outputDir;

  public static void main(String[] args) {
    CliHelper cliHelper = new CliHelper(MethodHandles.lookup().lookupClass())
        .addOption("o", "output-dir", "directory to output to", true, "o");
    try {
      if (!cliHelper.parse(args)) {
        System.exit(1);
      }
      Path outputDir = cliHelper.getValidDirectory("o", true);
      SyntheticBatchTest piplelineTest = new SyntheticBatchTest(outputDir);
      piplelineTest.execute();
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  private void execute() throws Exception {

    makeReport("cftr.ref_ref", new String[]{
        "DPYD/novariant.vcf",
        "UGT1A1/s1s1.vcf",
        "TPMT/s1s1.vcf",
        "CYP3A5/s1s7.vcf",
        "CFTR/refref.vcf",
        "cyp2c19/s2s2.vcf",
        "CYP2C9/s2s3.vcf",
        "SLCO1B1/s5s15.vcf",
        "VKORC1/-1639A-1639A.vcf",
        "cyp4f2/s1s1.vcf",
        "IFNL3/rs12979860CC.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("example", new String[]{
        "DPYD/novariant.vcf",
        "UGT1A1/s1s1.vcf",
        "TPMT/s1s1.vcf",
        "CYP3A5/s1s1.vcf",
        "CFTR/refref.vcf",
        "cyp2c19/s1s1.vcf",
        "CYP2C9/s1s1.vcf",
        "SLCO1B1/s1as1a.vcf",
        "VKORC1/-1639A-1639A.vcf",
        "cyp4f2/s1s1.vcf",
        "IFNL3/rs12979860CC.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("example.CYP2C19_noCYP2D6", new String[]{
        "cyp2c19/s1s1.vcf"
    }, null);



    makeReport("slco1b1.17.21", new String[]{
        "DPYD/novariant.vcf",
        "UGT1A1/s1s1.vcf",
        "TPMT/s1s1.vcf",
        "CYP3A5/s1s7.vcf",
        "CFTR/refref.vcf",
        "cyp2c19/s2s2.vcf",
        "CYP2C9/s2s3.vcf",
        "SLCO1B1/s17s21.vcf",
        "VKORC1/-1639A-1639A.vcf",
        "cyp4f2/s1s1.vcf",
        "IFNL3/rs12979860CC.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("slco1b1.5.15", new String[]{
        "DPYD/novariant.vcf",
        "UGT1A1/s1s1.vcf",
        "TPMT/s1s1.vcf",
        "CYP3A5/s1s7.vcf",
        "CFTR/refref.vcf",
        "cyp2c19/s2s2.vcf",
        "CYP2C9/s2s3.vcf",
        "SLCO1B1/s5s15.vcf",
        "VKORC1/-1639A-1639A.vcf",
        "cyp4f2/s1s1.vcf",
        "IFNL3/rs12979860CC.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("slco1b1.missing", new String[]{
        "DPYD/novariant.vcf",
        "UGT1A1/s1s1.vcf",
        "TPMT/s1s1.vcf",
        "CYP3A5/s1s7.vcf",
        "CFTR/refref.vcf",
        "cyp2c19/s2s2.vcf",
        "CYP2C9/s2s3.vcf",
        "VKORC1/-1639A-1639A.vcf",
        "cyp4f2/s1s1.vcf",
        "IFNL3/rs12979860CC.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("slco1b1.multi", new String[]{
        "DPYD/novariant.vcf",
        "UGT1A1/s1s1.vcf",
        "TPMT/s1s1.vcf",
        "CYP3A5/s1s7.vcf",
        "CFTR/refref.vcf",
        "cyp2c19/s2s2.vcf",
        "CYP2C9/s2s3.vcf",
        "SLCO1B1/multi.vcf",
        "VKORC1/-1639A-1639A.vcf",
        "cyp4f2/s1s1.vcf",
        "IFNL3/rs12979860CC.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("slco1b1.no_match", new String[]{
        "SLCO1B1/multi.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("slco1b1.match_1", new String[]{
        "SLCO1B1/s5s15.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("cyp2c19.onlyRs12769205", new String[]{
        "DPYD/novariant.vcf",
        "UGT1A1/s1s1.vcf",
        "TPMT/s1s1.vcf",
        "CYP3A5/s1s7.vcf",
        "CFTR/refref.vcf",
        "cyp2c19/rs12769205only.vcf",
        "CYP2C9/s1s1.vcf",
        "SLCO1B1/s1as1a.vcf",
        "VKORC1/-1639A-1639A.vcf",
        "cyp4f2/s1s1.vcf",
        "IFNL3/rs12979860CC.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("cyp2c19.refRs12769205", new String[]{
        "DPYD/novariant.vcf",
        "UGT1A1/s1s1.vcf",
        "TPMT/s1s1.vcf",
        "CYP3A5/s1s7.vcf",
        "CFTR/refref.vcf",
        "cyp2c19/rs12769205ref.vcf",
        "CYP2C9/s1s1.vcf",
        "SLCO1B1/s1as1a.vcf",
        "VKORC1/-1639A-1639A.vcf",
        "cyp4f2/s1s1.vcf",
        "IFNL3/rs12979860CC.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("cyp2c19.s2s3", new String[]{
        "DPYD/novariant.vcf",
        "UGT1A1/s1s1.vcf",
        "TPMT/s1s1.vcf",
        "CYP3A5/s1s7.vcf",
        "CFTR/refref.vcf",
        "cyp2c19/s2s3.vcf",
        "CYP2C9/s1s1.vcf",
        "SLCO1B1/s1as1a.vcf",
        "VKORC1/-1639A-1639A.vcf",
        "cyp4f2/s1s1.vcf",
        "IFNL3/rs12979860CC.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("cyp2c19.s4s17missingS1", new String[]{
        "cyp2c19/s4s17missingS1.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("cyp2c19.s1s4s17", new String[]{
        "cyp2c19/s1s4s17.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("ugt1a1.phased.multi", new String[]{
        "DPYD/novariant.vcf",
        "UGT1A1/s1s60s80phased.vcf",
        "TPMT/s1s1.vcf",
        "CYP3A5/s1s7.vcf",
        "cyp2c19/s2s2.vcf",
        "CYP2C9/s2s3.vcf",
        "SLCO1B1/s5s15.vcf",
        "VKORC1/-1639A-1639A.vcf",
        "cyp4f2/s1s1.vcf",
        "IFNL3/rs12979860CC.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("ugt1a1.unphased.multi", new String[]{
        "DPYD/novariant.vcf",
        "UGT1A1/s1s60s80unphased.vcf",
        "TPMT/s1s1.vcf",
        "CYP3A5/s1s7.vcf",
        "cyp2c19/s2s2.vcf",
        "CYP2C9/s2s3.vcf",
        "SLCO1B1/s5s15.vcf",
        "VKORC1/-1639A-1639A.vcf",
        "cyp4f2/s1s1.vcf",
        "IFNL3/rs12979860CC.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("missing.genes", new String[]{
        "DPYD/novariant.vcf",
        "TPMT/s1s1.vcf",
        "CYP3A5/s1s7.vcf",
        "cyp2c19/s2s2.vcf",
        "SLCO1B1/s5s15.vcf",
        "VKORC1/-1639A-1639A.vcf",
        "cyp4f2/s1s1.vcf",
        "IFNL3/rs12979860CC.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("cyp2c19.nonnormal", new String[]{
        "DPYD/novariant.vcf",
        "TPMT/s1s1.vcf",
        "CYP3A5/s1s7.vcf",
        "cyp2c19/s1s35.vcf",
        "SLCO1B1/s5s15.vcf",
        "VKORC1/-1639A-1639A.vcf",
        "cyp4f2/s1s1.vcf",
        "IFNL3/rs12979860CC.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("tpmt.star1s", new String[]{
        "DPYD/novariant.vcf",
        "UGT1A1/s1s1.vcf",
        "TPMT/s1ss1ss3.vcf",
        "CYP3A5/s1s7.vcf",
        "CFTR/refref.vcf",
        "cyp2c19/s2s2.vcf",
        "CYP2C9/s2s3.vcf",
        "SLCO1B1/s5s15.vcf",
        "VKORC1/-1639A-1639A.vcf",
        "cyp4f2/s1s1.vcf",
        "IFNL3/rs12979860CC.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("tpmt.s1s1s", new String[]{
        "DPYD/novariant.vcf",
        "UGT1A1/s1s1.vcf",
        "TPMT/s1s1s.vcf",
        "CYP3A5/s1s7.vcf",
        "CFTR/refref.vcf",
        "cyp2c19/s2s2.vcf",
        "CYP2C9/s2s3.vcf",
        "SLCO1B1/s5s15.vcf",
        "VKORC1/-1639A-1639A.vcf",
        "cyp4f2/s1s1.vcf",
        "IFNL3/rs12979860CC.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("tpmt.hom1s_het3a", new String[]{
        "DPYD/novariant.vcf",
        "UGT1A1/s1s1.vcf",
        "TPMT/hom1s_het3a.vcf",
        "CYP3A5/s1s7.vcf",
        "CFTR/refref.vcf",
        "cyp2c19/s2s2.vcf",
        "CYP2C9/s2s3.vcf",
        "SLCO1B1/s5s15.vcf",
        "VKORC1/-1639A-1639A.vcf",
        "cyp4f2/s1s1.vcf",
        "IFNL3/rs12979860CC.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("tpmt.het3a", new String[]{
        "DPYD/novariant.vcf",
        "UGT1A1/s1s1.vcf",
        "TPMT/het3a.vcf",
        "CYP3A5/s1s7.vcf",
        "CFTR/refref.vcf",
        "cyp2c19/s2s2.vcf",
        "CYP2C9/s2s3.vcf",
        "SLCO1B1/s5s15.vcf",
        "VKORC1/-1639A-1639A.vcf",
        "cyp4f2/s1s1.vcf",
        "IFNL3/rs12979860CC.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("tpmt.s15offdata", new String[]{
        "DPYD/novariant.vcf",
        "UGT1A1/s1s1.vcf",
        "TPMT/s15offdata.vcf",
        "CYP3A5/s1s7.vcf",
        "CFTR/refref.vcf",
        "cyp2c19/s2s2.vcf",
        "CYP2C9/s2s3.vcf",
        "SLCO1B1/s5s15.vcf",
        "VKORC1/-1639A-1639A.vcf",
        "cyp4f2/s1s1.vcf",
        "IFNL3/rs12979860CC.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("tpmt.unexpected.allele", new String[]{
        "TPMT/unexpected.allele.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("tpmt.unexpected.allele.with.s9", new String[]{
        "TPMT/unexpected.allele.with.s9.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("cyp2c19.rs28399504missing", new String[]{
        "DPYD/novariant.vcf",
        "UGT1A1/s1s1.vcf",
        "TPMT/s1ss1ss3.vcf",
        "CYP3A5/s1s7.vcf",
        "CFTR/refref.vcf",
        "cyp2c19/s4bs17rs28399504missing.vcf",
        "CYP2C9/s2s3.vcf",
        "SLCO1B1/s5s15.vcf",
        "VKORC1/-1639A-1639A.vcf",
        "cyp4f2/s1s1.vcf",
        "IFNL3/rs12979860CC.vcf"
    }, sf_outsideCYP2D6File);


    makeReport("dpyd.stars12b", new String[]{
        "DPYD/s1s2b.vcf",
        "UGT1A1/s1s1.vcf",
        "TPMT/s1ss1ss3.vcf",
        "CYP3A5/s1s7.vcf",
        "CFTR/refref.vcf",
        "cyp2c19/s2s2.vcf",
        "CYP2C9/s2s3.vcf",
        "SLCO1B1/s5s15.vcf",
        "VKORC1/-1639A-1639A.vcf",
        "cyp4f2/s1s1.vcf",
        "IFNL3/rs12979860CC.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("cyp2c19.rxPossible", new String[]{
        "cyp2c19/s1s17.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("cyp2c19.s2s2", new String[]{
        "cyp2c19/s2s2.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("cyp2c19.s1s2", new String[]{
        "cyp2c19/s1s2.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("cyp2c19.s1s2rs3758581missing", new String[]{
        "cyp2c19/s1s2rs3758581missing.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("cyp2c19.s1s2rs3758581missingRs58973490het", new String[]{
        "cyp2c19/s1s2rs3758581missingRs58973490het.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("cyp2c19.s1s2rs58973490het", new String[]{
        "cyp2c19/s1s2rs58973490het.vcf"
    }, sf_outsideCYP2D6File);


    makeReport("ugt1a1.s1s1", new String[]{
        "UGT1A1/s1s1.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("ugt1a1.s1s28s60s80unphased", new String[]{
        "UGT1A1/s1s28s60s80unphased.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("ugt1a1.s1s60s80unphased", new String[]{
        "UGT1A1/s1s60s80unphased.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("ugt1a1.s1s60s80phased", new String[]{
        "UGT1A1/s1s60s80phased.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("ugt1a1.s6s6", new String[]{
        "UGT1A1/s6s6.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("ugt1a1.s28s37", new String[]{
        "UGT1A1/s28s37.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("ugt1a1.s28s80unphased", new String[]{
        "UGT1A1/s28s80unphased.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("ugt1a1.s28s80phased", new String[]{
        "UGT1A1/s28s80phased.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("ugt1a1.phased.s28s80s6s60", new String[]{
        "UGT1A1/s28s80s6s60phased.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("ugt1a1.unphased.s28s80s6s60", new String[]{
        "UGT1A1/s28s80s6s60unphased.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("ugt1a1.s6s60s80s28missingphased", new String[]{
        "UGT1A1/s6s60s80s28missingphased.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("ugt1a1.s6s60s80s28missingunphased", new String[]{
        "UGT1A1/s6s60s80s28missingunphased.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("ugt1a1.s80s28missing", new String[]{
        "UGT1A1/s80s28missing.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("ugt1a1.s28het", new String[]{
        "UGT1A1/s28het.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("ugt1a1.na12717", new String[]{
        "UGT1A1/NA12717_UGT1A1.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("ugt1a1.na18868", new String[]{
        "UGT1A1/NA18868_UGT1A1.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("ugt1a1.na19785", new String[]{
        "UGT1A1/NA19785_UGT1A1.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("ugt1a1.s28s28unphaseds60s80miss", new String[]{
        "UGT1A1/s28s28unphaseds60s80miss.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("ugt1a1.s27s28unphaseds80s60missing", new String[]{
        "UGT1A1/s27s28unphaseds80s60missing.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("ugt1a1.s28hets60homounphaseds80missing", new String[]{
        "UGT1A1/s28hets60homounphaseds80missing.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("ugt1a1.HG00436", new String[]{
        "UGT1A1/HG00436.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("ugt1a1.s1s80s27s60s28missingphased", new String[]{
        "UGT1A1/s1s80s27s60s28missingphased.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("ugt1a1.s1s60s80s6phased", new String[]{
        "UGT1A1/s1s60s80s6phased.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("ugt1a1.s1s60s80s28s6phased", new String[]{
        "UGT1A1/s1s60s80s28s6phased.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("ugt1a1.s1s37s80s60phased", new String[]{
        "UGT1A1/s1s37s80s60phased.vcf"
    }, sf_outsideCYP2D6File);


    makeReport("cyp2c19.rs12248560missing", new String[]{
        "cyp2c19/s1s1rs12248560missing.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("cyp3a5.s1s1rs776746missing", new String[]{
        "cyp3a5/s1s1rs776746missing.vcf"
    }, sf_outsideCYP2D6File);
    makeReport("cyp3a5.s1s3rs776746rs55965422het", new String[]{
        "cyp3a5/s1s3rs776746rs55965422het.vcf"
    }, sf_outsideCYP2D6File);
    makeReport("cyp3a5.s1s3rs776746rs55965422rs28383479het", new String[]{
        "cyp3a5/s1s3rs776746rs55965422rs28383479het.vcf"
    }, sf_outsideCYP2D6File);
    makeReport("cyp3a5.s3s3rs55965422het", new String[]{
        "cyp3a5/s3s3rs55965422het.vcf"
    }, sf_outsideCYP2D6File);
    makeReport("cyp3a5.s3s5-homozygous", new String[]{
        "cyp3a5/s3s5-homozygous.vcf"
    }, sf_outsideCYP2D6File);
    makeReport("cyp3a5.s1s3rs776746rs28383479het", new String[]{
        "cyp3a5/s1s3rs776746rs28383479het.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("cftr.G1244Eref", new String[]{
        "cftr/G1244Eref.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("cftr.ref_vwithrecom", new String[]{
        "cftr/ref_vwithrecom.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("cftr.v_vwithrecom", new String[]{
        "cftr/v_vwithrecom.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("cftr.v_vwithrecom1", new String[]{
        "cftr/v_vwithrecom1.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("cyp2c9.s2s24", new String[]{
        "cyp2c9/s2s24.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("cftr.refI507missing", new String[]{
        "cftr/refI507missing.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("cyp2c9.s2s24only", new String[]{
        "cyp2c9/s2s24only.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("cyp2c9.s1s61", new String[]{
        "cyp2c9/s1s61.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("cyp2c9.s3s18m94981296", new String[]{
        "cyp2c9/s3s18m94981296.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("cyp2c19.rs12769205call", new String[]{
        "cyp2c19/rs12769205call.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("cyp2c19.s4s17het", new String[]{
        "cyp2c19/s4s17het.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("dpyd.c557c703", new String[]{
        "DPYD/c557c703.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("dpyd.c1129-5923c2846", new String[]{
        "DPYD/c1129-5923c2846.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("dpyd.c1156het", new String[]{
        "DPYD/c1156het.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("dpyd.c1679c1156", new String[]{
        "DPYD/c1679c1156.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("dpyd.c2846het", new String[]{
        "DPYD/c2846het.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("dpyd.novariant", new String[]{
        "DPYD/novariant.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("nudt15.refref", new String[]{
        "NUDT15/refref.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("nudt15.s3ref", new String[]{
        "NUDT15/s3ref.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("g6pd.homB", new String[]{
        "DPYD/novariant.vcf",
        "UGT1A1/s1s1.vcf",
        "TPMT/s1s1.vcf",
        "CYP3A5/s1s1.vcf",
        "CFTR/refref.vcf",
        "cyp2c19/s1s1.vcf",
        "CYP2C9/s1s1.vcf",
        "SLCO1B1/s1as1a.vcf",
        "VKORC1/-1639A-1639A.vcf",
        "cyp4f2/s1s1.vcf",
        "IFNL3/rs12979860CC.vcf"
    }, sf_outsideCYP2D6G6PDFile);

    makeReport("cyp2b6.s1s1", new String[]{
        "CYP2B6/s1s1.vcf"
    }, null);

    makeReport("cyp2b6.s1s34", new String[]{
        "CYP2B6/s1s34.vcf"
    }, null);

    makeReportWithOutputString("cyp2d6.s1unknown", new String[]{
        "cyp2c9/s1s61.vcf"
    }, "CYP2D6\t*1/*XXX\n");

    makeReportWithOutputString("mtrnr1.increased", new String[]{
        "cyp2c9/s1s61.vcf"
    }, "MT-RNR1\t1555A>G\n");

    sf_logger.info("Wrote reports to {}", f_outputDir);
  }


  private SyntheticBatchTest(Path outputDir) throws IOException {
    f_outputDir = outputDir;
    f_pharmcat = new PharmCAT(outputDir, null).keepMatcherOutput();
    f_pharmcat
        .writeJson(true)
        .writePhenotyperJson(true);
  }

  private void makeReport(String key, String[] testVcfs, Path outsideCallPath) throws Exception {
    Path sampleDir = f_outputDir.resolve(key);
    if (!sampleDir.toFile().exists()) {
      if (!sampleDir.toFile().mkdirs()) {
        throw new RuntimeException("Output directory could not be created " + sampleDir.toAbsolutePath());
      }
    }

    Path sampleVcf = writeVcf(sampleDir.resolve(key + ".vcf"), testVcfs);
    f_pharmcat.setOutputDir(sampleDir);
    f_pharmcat.execute(sampleVcf, outsideCallPath, null);
  }

  private void makeReportWithOutputString(String key, String[] testVcfs, String outsideCalls) throws Exception {
    Path outsideCallPath = Files.createTempFile("outsideCall", ".tsv");
    try (FileWriter fw = new FileWriter(outsideCallPath.toFile())) {
      fw.write(outsideCalls);
    }
    makeReport(key, testVcfs, outsideCallPath);
    outsideCallPath.toFile().deleteOnExit();
  }

  private Path writeVcf(Path outputVcf, String[] filesToInclude) {
    try (FileWriter writer = new FileWriter(outputVcf.toFile())) {
      writer.write(VcfTestUtils.writeVcf(filesToInclude));
    } catch (IOException e) {
      e.printStackTrace();
    }

    return outputVcf;
  }
}

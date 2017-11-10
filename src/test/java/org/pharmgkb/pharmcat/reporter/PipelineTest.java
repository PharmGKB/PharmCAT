package org.pharmgkb.pharmcat.reporter;

import java.io.FileWriter;
import java.io.IOException;
import java.lang.invoke.MethodHandles;
import java.nio.file.Path;
import java.util.Date;
import java.util.LinkedHashMap;
import java.util.Map;
import org.pharmgkb.common.io.util.CliHelper;
import org.pharmgkb.common.util.PathUtils;
import org.pharmgkb.pharmcat.PharmCAT;
import org.pharmgkb.pharmcat.VcfTestUtils;


/**
 * This test generates a collection of synthetic sample VCF files from our test resources and then runs the matcher and
 * reporter on those samples, writing the output to a desired directory.
 *
 * This is especially useful to see how different allele calls will be called and displayed in the final output report.
 *
 * @author Ryan Whaley
 */
public class PipelineTest {
  private static final Path sf_astrolabe = PathUtils.getPathToResource("org/pharmgkb/pharmcat/reporter/test.astrolabe.tsv");
  private static final Map<String,String[]> sf_testVcfs = new LinkedHashMap<>();
  static {

    sf_testVcfs.put("test.cftr.reg_inc", new String[]{
        "DPYD/s1s1.vcf",
        "UGT1A1/s1s1.vcf",
        "TPMT/s1s1.vcf",
        "CYP3A5/s1s7.vcf",
        "CFTR/G542XF508del.vcf",
        "CYP2C19/s2s2.vcf",
        "CYP2C9/s2s3.vcf",
        "SLCO1B1/s5s15.vcf",
        "VKORC1/-1639A-1639A.vcf",
        "cyp4f2/s1s1.vcf",
        "IFNL3/rs12979860CC.vcf"
    });

    sf_testVcfs.put("test.cftr.ref_inc", new String[]{
        "DPYD/s1s1.vcf",
        "UGT1A1/s1s1.vcf",
        "TPMT/s1s1.vcf",
        "CYP3A5/s1s7.vcf",
        "CFTR/refG542X.vcf",
        "CYP2C19/s2s2.vcf",
        "CYP2C9/s2s3.vcf",
        "SLCO1B1/s5s15.vcf",
        "VKORC1/-1639A-1639A.vcf",
        "cyp4f2/s1s1.vcf",
        "IFNL3/rs12979860CC.vcf"
    });

    sf_testVcfs.put("test.cftr.inc_inc", new String[]{
        "DPYD/s1s1.vcf",
        "UGT1A1/s1s1.vcf",
        "TPMT/s1s1.vcf",
        "CYP3A5/s1s7.vcf",
        "CFTR/G542XG542X.vcf",
        "CYP2C19/s2s2.vcf",
        "CYP2C9/s2s3.vcf",
        "SLCO1B1/s5s15.vcf",
        "VKORC1/-1639A-1639A.vcf",
        "cyp4f2/s1s1.vcf",
        "IFNL3/rs12979860CC.vcf"
    });

    sf_testVcfs.put("test.cftr.ref_ref", new String[]{
        "DPYD/s1s1.vcf",
        "UGT1A1/s1s1.vcf",
        "TPMT/s1s1.vcf",
        "CYP3A5/s1s7.vcf",
        "CFTR/refref.vcf",
        "CYP2C19/s2s2.vcf",
        "CYP2C9/s2s3.vcf",
        "SLCO1B1/s5s15.vcf",
        "VKORC1/-1639A-1639A.vcf",
        "cyp4f2/s1s1.vcf",
        "IFNL3/rs12979860CC.vcf"
    });

    sf_testVcfs.put("test.cftr.reg_reg", new String[]{
        "DPYD/s1s1.vcf",
        "UGT1A1/s1s1.vcf",
        "TPMT/s1s1.vcf",
        "CYP3A5/s1s7.vcf",
        "CFTR/F508delF508del.vcf",
        "CYP2C19/s2s2.vcf",
        "CYP2C9/s2s3.vcf",
        "SLCO1B1/s5s15.vcf",
        "VKORC1/-1639A-1639A.vcf",
        "cyp4f2/s1s1.vcf",
        "IFNL3/rs12979860CC.vcf"
    });

    sf_testVcfs.put("test.cftr.ref_reg", new String[]{
        "DPYD/s1s1.vcf",
        "UGT1A1/s1s1.vcf",
        "TPMT/s1s1.vcf",
        "CYP3A5/s1s7.vcf",
        "CFTR/refF508del.vcf",
        "CYP2C19/s2s2.vcf",
        "CYP2C9/s2s3.vcf",
        "SLCO1B1/s5s15.vcf",
        "VKORC1/-1639A-1639A.vcf",
        "cyp4f2/s1s1.vcf",
        "IFNL3/rs12979860CC.vcf"
    });



    sf_testVcfs.put("test.slco1b1.17.21", new String[]{
        "DPYD/s1s1.vcf",
        "UGT1A1/s1s1.vcf",
        "TPMT/s1s1.vcf",
        "CYP3A5/s1s7.vcf",
        "CFTR/refref.vcf",
        "CYP2C19/s2s2.vcf",
        "CYP2C9/s2s3.vcf",
        "SLCO1B1/s17s21.vcf",
        "VKORC1/-1639A-1639A.vcf",
        "cyp4f2/s1s1.vcf",
        "IFNL3/rs12979860CC.vcf"
    });

    sf_testVcfs.put("test.slco1b1.5.15", new String[]{
        "DPYD/s1s1.vcf",
        "UGT1A1/s1s1.vcf",
        "TPMT/s1s1.vcf",
        "CYP3A5/s1s7.vcf",
        "CFTR/refref.vcf",
        "CYP2C19/s2s2.vcf",
        "CYP2C9/s2s3.vcf",
        "SLCO1B1/s5s15.vcf",
        "VKORC1/-1639A-1639A.vcf",
        "cyp4f2/s1s1.vcf",
        "IFNL3/rs12979860CC.vcf"
    });

    sf_testVcfs.put("test.slco1b1.missing", new String[]{
        "DPYD/s1s1.vcf",
        "UGT1A1/s1s1.vcf",
        "TPMT/s1s1.vcf",
        "CYP3A5/s1s7.vcf",
        "CFTR/refref.vcf",
        "CYP2C19/s2s2.vcf",
        "CYP2C9/s2s3.vcf",
        "VKORC1/-1639A-1639A.vcf",
        "cyp4f2/s1s1.vcf",
        "IFNL3/rs12979860CC.vcf"
    });

    sf_testVcfs.put("test.slco1b1.multi", new String[]{
        "DPYD/s1s1.vcf",
        "UGT1A1/s1s1.vcf",
        "TPMT/s1s1.vcf",
        "CYP3A5/s1s7.vcf",
        "CFTR/refref.vcf",
        "CYP2C19/s2s2.vcf",
        "CYP2C9/s2s3.vcf",
        "SLCO1B1/multi.vcf",
        "VKORC1/-1639A-1639A.vcf",
        "cyp4f2/s1s1.vcf",
        "IFNL3/rs12979860CC.vcf"
    });

    sf_testVcfs.put("test.slco1b1.no_match", new String[]{
        "SLCO1B1/multi.vcf"
    });

    sf_testVcfs.put("test.slco1b1.match_1", new String[]{
        "SLCO1B1/s5s15.vcf"
    });

    sf_testVcfs.put("test.cyp2c19.onlyRs12769205", new String[]{
        "DPYD/s1s1.vcf",
        "UGT1A1/s1s1.vcf",
        "TPMT/s1s1.vcf",
        "CYP3A5/s1s7.vcf",
        "CFTR/refref.vcf",
        "CYP2C19/rs12769205only.vcf",
        "CYP2C9/s1s1.vcf",
        "SLCO1B1/s1as1a.vcf",
        "VKORC1/-1639A-1639A.vcf",
        "cyp4f2/s1s1.vcf",
        "IFNL3/rs12979860CC.vcf"
    });

    sf_testVcfs.put("test.cyp2c19.refRs12769205", new String[]{
        "DPYD/s1s1.vcf",
        "UGT1A1/s1s1.vcf",
        "TPMT/s1s1.vcf",
        "CYP3A5/s1s7.vcf",
        "CFTR/refref.vcf",
        "CYP2C19/rs12769205ref.vcf",
        "CYP2C9/s1s1.vcf",
        "SLCO1B1/s1as1a.vcf",
        "VKORC1/-1639A-1639A.vcf",
        "cyp4f2/s1s1.vcf",
        "IFNL3/rs12979860CC.vcf"
    });

    sf_testVcfs.put("test.cyp2c19.s2s3", new String[]{
        "DPYD/s1s1.vcf",
        "UGT1A1/s1s1.vcf",
        "TPMT/s1s1.vcf",
        "CYP3A5/s1s7.vcf",
        "CFTR/refref.vcf",
        "CYP2C19/s2s3.vcf",
        "CYP2C9/s1s1.vcf",
        "SLCO1B1/s1as1a.vcf",
        "VKORC1/-1639A-1639A.vcf",
        "cyp4f2/s1s1.vcf",
        "IFNL3/rs12979860CC.vcf"
    });

    sf_testVcfs.put("test.ugt1a1.phased.multi", new String[]{
        "DPYD/s1s1.vcf",
        "UGT1A1/s1s60s80phased.vcf",
        "TPMT/s1s1.vcf",
        "CYP3A5/s1s7.vcf",
        "CFTR/F508delF508del.vcf",
        "CYP2C19/s2s2.vcf",
        "CYP2C9/s2s3.vcf",
        "SLCO1B1/s5s15.vcf",
        "VKORC1/-1639A-1639A.vcf",
        "cyp4f2/s1s1.vcf",
        "IFNL3/rs12979860CC.vcf"
    });

    sf_testVcfs.put("test.ugt1a1.unphased.multi", new String[]{
        "DPYD/s1s1.vcf",
        "UGT1A1/s1s60s80.vcf",
        "TPMT/s1s1.vcf",
        "CYP3A5/s1s7.vcf",
        "CFTR/F508delF508del.vcf",
        "CYP2C19/s2s2.vcf",
        "CYP2C9/s2s3.vcf",
        "SLCO1B1/s5s15.vcf",
        "VKORC1/-1639A-1639A.vcf",
        "cyp4f2/s1s1.vcf",
        "IFNL3/rs12979860CC.vcf"
    });

    sf_testVcfs.put("test.missing.genes", new String[]{
        "DPYD/s1s1.vcf",
        "TPMT/s1s1.vcf",
        "CYP3A5/s1s7.vcf",
        "CFTR/F508delF508del.vcf",
        "CYP2C19/s2s2.vcf",
        "SLCO1B1/s5s15.vcf",
        "VKORC1/-1639A-1639A.vcf",
        "cyp4f2/s1s1.vcf",
        "IFNL3/rs12979860CC.vcf"
    });

    sf_testVcfs.put("test.cyp2c19.nonnormal", new String[]{
        "DPYD/s1s1.vcf",
        "TPMT/s1s1.vcf",
        "CYP3A5/s1s7.vcf",
        "CFTR/F508delF508del.vcf",
        "CYP2C19/s1s35.vcf",
        "SLCO1B1/s5s15.vcf",
        "VKORC1/-1639A-1639A.vcf",
        "cyp4f2/s1s1.vcf",
        "IFNL3/rs12979860CC.vcf"
    });

    sf_testVcfs.put("test.tpmt.star1s", new String[]{
        "DPYD/s1s1.vcf",
        "UGT1A1/s1s1.vcf",
        "TPMT/s1ss1ss3.vcf",
        "CYP3A5/s1s7.vcf",
        "CFTR/G542XF508del.vcf",
        "CYP2C19/s2s2.vcf",
        "CYP2C9/s2s3.vcf",
        "SLCO1B1/s5s15.vcf",
        "VKORC1/-1639A-1639A.vcf",
        "cyp4f2/s1s1.vcf",
        "IFNL3/rs12979860CC.vcf"
    });

    sf_testVcfs.put("test.tpmt.s1s1s", new String[]{
        "DPYD/s1s1.vcf",
        "UGT1A1/s1s1.vcf",
        "TPMT/s1s1s.vcf",
        "CYP3A5/s1s7.vcf",
        "CFTR/G542XF508del.vcf",
        "CYP2C19/s2s2.vcf",
        "CYP2C9/s2s3.vcf",
        "SLCO1B1/s5s15.vcf",
        "VKORC1/-1639A-1639A.vcf",
        "cyp4f2/s1s1.vcf",
        "IFNL3/rs12979860CC.vcf"
    });

    sf_testVcfs.put("test.tpmt.hom1s_het3a", new String[]{
        "DPYD/s1s1.vcf",
        "UGT1A1/s1s1.vcf",
        "TPMT/hom1s_het3a.vcf",
        "CYP3A5/s1s7.vcf",
        "CFTR/G542XF508del.vcf",
        "CYP2C19/s2s2.vcf",
        "CYP2C9/s2s3.vcf",
        "SLCO1B1/s5s15.vcf",
        "VKORC1/-1639A-1639A.vcf",
        "cyp4f2/s1s1.vcf",
        "IFNL3/rs12979860CC.vcf"
    });

    sf_testVcfs.put("test.tpmt.het3a", new String[]{
        "DPYD/s1s1.vcf",
        "UGT1A1/s1s1.vcf",
        "TPMT/het3a.vcf",
        "CYP3A5/s1s7.vcf",
        "CFTR/G542XF508del.vcf",
        "CYP2C19/s2s2.vcf",
        "CYP2C9/s2s3.vcf",
        "SLCO1B1/s5s15.vcf",
        "VKORC1/-1639A-1639A.vcf",
        "cyp4f2/s1s1.vcf",
        "IFNL3/rs12979860CC.vcf"
    });

    sf_testVcfs.put("test.cyp2c19.rs28399504missing", new String[]{
        "DPYD/s1s1.vcf",
        "UGT1A1/s1s1.vcf",
        "TPMT/s1ss1ss3.vcf",
        "CYP3A5/s1s7.vcf",
        "CFTR/G542XF508del.vcf",
        "CYP2C19/s4bs17rs28399504missing.vcf",
        "CYP2C9/s2s3.vcf",
        "SLCO1B1/s5s15.vcf",
        "VKORC1/-1639A-1639A.vcf",
        "cyp4f2/s1s1.vcf",
        "IFNL3/rs12979860CC.vcf"
    });


    sf_testVcfs.put("test.dpyd.stars12b", new String[]{
        "DPYD/s1s2b.vcf",
        "UGT1A1/s1s1.vcf",
        "TPMT/s1ss1ss3.vcf",
        "CYP3A5/s1s7.vcf",
        "CFTR/G542XF508del.vcf",
        "CYP2C19/s2s2.vcf",
        "CYP2C9/s2s3.vcf",
        "SLCO1B1/s5s15.vcf",
        "VKORC1/-1639A-1639A.vcf",
        "cyp4f2/s1s1.vcf",
        "IFNL3/rs12979860CC.vcf"
    });

    sf_testVcfs.put("test.cyp2c19.rxPossible", new String[]{
        "CYP2C19/s1s17.vcf"
    });


    sf_testVcfs.put("test.ugt1a1.s1s1", new String[]{
      "UGT1A1/s1s1.vcf"
    });

    sf_testVcfs.put("test.ugt1a1.s1s28s60s80", new String[]{
      "UGT1A1/s1s28s60s80.vcf"
    });

    sf_testVcfs.put("test.ugt1a1.s1s60s80", new String[]{
      "UGT1A1/s1s60s80.vcf"
    });

    sf_testVcfs.put("test.ugt1a1.s1s60s80phased", new String[]{
      "UGT1A1/s1s60s80phased.vcf"
    });

    sf_testVcfs.put("test.ugt1a1.s6s6", new String[]{
      "UGT1A1/s6s6.vcf"
    });

    sf_testVcfs.put("test.ugt1a1.s28s37", new String[]{
      "UGT1A1/s28s37.vcf"
    });

    sf_testVcfs.put("test.ugt1a1.s28s80", new String[]{
      "UGT1A1/s28s80.vcf"
    });

    sf_testVcfs.put("test.ugt1a1.phased.s28s80s6s60", new String[]{
        "UGT1A1/s28s80s6s60phased.vcf"
    });

    sf_testVcfs.put("test.ugt1a1.unphased.s28s80s6s60", new String[]{
        "UGT1A1/s28s80s6s60unphased.vcf"
    });

    sf_testVcfs.put("test.cyp2c19.rs12248560missing", new String[]{
        "cyp2c19/s1s1rs12248560missing.vcf"
    });

    sf_testVcfs.put("test.cyp3a5.s1s1rs776746missing", new String[]{
        "cyp3a5/s1s1rs776746missing.vcf"
    });

    sf_testVcfs.put("test.cftr.G1244Eref", new String[]{
        "cftr/G1244Eref.vcf"
    });

    sf_testVcfs.put("test.cftr.G1244EF508del", new String[]{
        "cftr/G1244EF508del.vcf"
    });

    sf_testVcfs.put("test.cftr.G551DG542X", new String[]{
        "cftr/G551DG542X.vcf"
    });

    sf_testVcfs.put("test.cyp2c9.s2s24", new String[]{
        "cyp2c9/s2s24.vcf"
    });

    sf_testVcfs.put("test.cyp2c19.rs12769205call", new String[]{
        "cyp2c19/rs12769205call.vcf"
    });

    sf_testVcfs.put("test.cyp2c19.s4s17het", new String[]{
        "cyp2c19/s4s17het.vcf"
    });
  }

  private PharmCAT m_pharmcat;
  private Path m_outputDir;

  public static void main(String[] args) {
    CliHelper cliHelper = new CliHelper(MethodHandles.lookup().lookupClass())
        .addOption("o", "output-dir", "directory to output to", true, "o")
        .addOption("g", "guideline-dir", "directory of guideline annotations (JSON files)", false, "n");

    try {
      if (!cliHelper.parse(args)) {
        System.exit(1);
      }

      Path outputDir = cliHelper.getValidDirectory("o", true);
      Path guidelineDir = null;
      if (cliHelper.hasOption("g")) {
        guidelineDir = cliHelper.getValidDirectory("g", false);
      }

      PipelineTest piplelineTest = new PipelineTest(outputDir, guidelineDir);
      piplelineTest.execute();
    } catch (Exception e) {
      e.printStackTrace();
    }
  }


  private PipelineTest(Path outputDir, Path guidelineDir) throws IOException {
    m_outputDir = outputDir;
    m_pharmcat = new PharmCAT(outputDir, null, guidelineDir).keepMatcherOutput();
  }

  private void execute() throws Exception {
    System.out.println("Run time: " + new Date());

    for (String key : sf_testVcfs.keySet()) {

      Path sampleVcf = writeVcf(m_outputDir.resolve(key+".vcf"), sf_testVcfs.get(key));
      m_pharmcat.execute(sampleVcf, sf_astrolabe, null);
      System.out.println("Generated "+key);

    }
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

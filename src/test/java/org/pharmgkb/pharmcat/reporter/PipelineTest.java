package org.pharmgkb.pharmcat.reporter;

import java.io.FileWriter;
import java.io.IOException;
import java.lang.invoke.MethodHandles;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Collection;
import java.util.Date;
import com.google.common.collect.LinkedHashMultimap;
import com.google.common.collect.Multimap;
import org.pharmgkb.common.io.util.CliHelper;
import org.pharmgkb.common.util.PathUtils;
import org.pharmgkb.pharmcat.PharmCAT;


/**
 * This test generates a collection of synthetic sample VCF files from our test resources and then runs the matcher and
 * reporter on those samples, writing the output to a desired directory.
 *
 * This is especially useful to see how different allele calls will be called and displayed in the final output report.
 *
 * @author Ryan Whaley
 */
public class PipelineTest {
  private static final String sf_headerFile = "org/pharmgkb/pharmcat/haplotype/DPYD/s1s1.vcf";
  private static final Path sf_astrolabe = PathUtils.getPathToResource("org/pharmgkb/pharmcat/reporter/test.astrolabe.tsv");
  private static final Multimap<String,String> sf_testVcfs = LinkedHashMultimap.create();
  static {
    String key;

    key = "test.cftr.reg_inc";
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/DPYD/s1s1.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/UGT1A1/s1s1.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/TPMT/s1s1.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/CYP3A5/s1s7.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/CFTR/G542XF508del.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/CYP2C19/s2s2.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/CYP2C9/s2s3.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/SLCO1B1/s5s15.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/VKORC1/-1639A-1639A.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/cyp4f2/s1s1.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/IFNL3/rs12979860CC.vcf");

    key = "test.cftr.ref_inc";
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/DPYD/s1s1.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/UGT1A1/s1s1.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/TPMT/s1s1.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/CYP3A5/s1s7.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/CFTR/refG542X.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/CYP2C19/s2s2.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/CYP2C9/s2s3.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/SLCO1B1/s5s15.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/VKORC1/-1639A-1639A.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/cyp4f2/s1s1.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/IFNL3/rs12979860CC.vcf");

    key = "test.cftr.inc_inc";
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/DPYD/s1s1.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/UGT1A1/s1s1.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/TPMT/s1s1.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/CYP3A5/s1s7.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/CFTR/G542XG542X.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/CYP2C19/s2s2.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/CYP2C9/s2s3.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/SLCO1B1/s5s15.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/VKORC1/-1639A-1639A.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/cyp4f2/s1s1.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/IFNL3/rs12979860CC.vcf");

    key = "test.cftr.ref_ref";
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/DPYD/s1s1.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/UGT1A1/s1s1.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/TPMT/s1s1.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/CYP3A5/s1s7.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/CFTR/refref.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/CYP2C19/s2s2.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/CYP2C9/s2s3.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/SLCO1B1/s5s15.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/VKORC1/-1639A-1639A.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/cyp4f2/s1s1.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/IFNL3/rs12979860CC.vcf");

    key = "test.cftr.reg_reg";
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/DPYD/s1s1.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/UGT1A1/s1s1.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/TPMT/s1s1.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/CYP3A5/s1s7.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/CFTR/F508delF508del.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/CYP2C19/s2s2.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/CYP2C9/s2s3.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/SLCO1B1/s5s15.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/VKORC1/-1639A-1639A.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/cyp4f2/s1s1.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/IFNL3/rs12979860CC.vcf");

    key = "test.cftr.ref_reg";
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/DPYD/s1s1.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/UGT1A1/s1s1.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/TPMT/s1s1.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/CYP3A5/s1s7.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/CFTR/refF508del.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/CYP2C19/s2s2.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/CYP2C9/s2s3.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/SLCO1B1/s5s15.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/VKORC1/-1639A-1639A.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/cyp4f2/s1s1.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/IFNL3/rs12979860CC.vcf");



    key = "test.slco1b1.17.21";
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/DPYD/s1s1.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/UGT1A1/s1s1.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/TPMT/s1s1.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/CYP3A5/s1s7.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/CFTR/refref.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/CYP2C19/s2s2.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/CYP2C9/s2s3.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/SLCO1B1/s17s21.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/VKORC1/-1639A-1639A.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/cyp4f2/s1s1.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/IFNL3/rs12979860CC.vcf");

    key = "test.slco1b1.5.15";
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/DPYD/s1s1.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/UGT1A1/s1s1.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/TPMT/s1s1.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/CYP3A5/s1s7.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/CFTR/refref.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/CYP2C19/s2s2.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/CYP2C9/s2s3.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/SLCO1B1/s5s15.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/VKORC1/-1639A-1639A.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/cyp4f2/s1s1.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/IFNL3/rs12979860CC.vcf");

    key = "test.slco1b1.missing";
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/DPYD/s1s1.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/UGT1A1/s1s1.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/TPMT/s1s1.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/CYP3A5/s1s7.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/CFTR/refref.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/CYP2C19/s2s2.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/CYP2C9/s2s3.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/VKORC1/-1639A-1639A.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/cyp4f2/s1s1.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/IFNL3/rs12979860CC.vcf");

    key = "test.slco1b1.multi";
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/DPYD/s1s1.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/UGT1A1/s1s1.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/TPMT/s1s1.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/CYP3A5/s1s7.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/CFTR/refref.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/CYP2C19/s2s2.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/CYP2C9/s2s3.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/SLCO1B1/multi.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/VKORC1/-1639A-1639A.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/cyp4f2/s1s1.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/IFNL3/rs12979860CC.vcf");

    key = "test.cyp2c19.onlyRs12769205";
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/DPYD/s1s1.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/UGT1A1/s1s1.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/TPMT/s1s1.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/CYP3A5/s1s7.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/CFTR/refref.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/CYP2C19/rs12769205only.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/CYP2C9/s1s1.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/SLCO1B1/s1as1a.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/VKORC1/-1639A-1639A.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/cyp4f2/s1s1.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/IFNL3/rs12979860CC.vcf");

    key = "test.cyp2c19.refRs12769205";
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/DPYD/s1s1.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/UGT1A1/s1s1.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/TPMT/s1s1.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/CYP3A5/s1s7.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/CFTR/refref.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/CYP2C19/rs12769205ref.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/CYP2C9/s1s1.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/SLCO1B1/s1as1a.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/VKORC1/-1639A-1639A.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/cyp4f2/s1s1.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/IFNL3/rs12979860CC.vcf");

    key = "test.ugt1a1.phased.multi";
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/DPYD/s1s1.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/UGT1A1/s1s60s80phased.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/TPMT/s1s1.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/CYP3A5/s1s7.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/CFTR/F508delF508del.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/CYP2C19/s2s2.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/CYP2C9/s2s3.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/SLCO1B1/s5s15.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/VKORC1/-1639A-1639A.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/cyp4f2/s1s1.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/IFNL3/rs12979860CC.vcf");

    key = "test.ugt1a1.unphased.multi";
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/DPYD/s1s1.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/UGT1A1/s1s60s80.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/TPMT/s1s1.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/CYP3A5/s1s7.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/CFTR/F508delF508del.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/CYP2C19/s2s2.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/CYP2C9/s2s3.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/SLCO1B1/s5s15.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/VKORC1/-1639A-1639A.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/cyp4f2/s1s1.vcf");
    sf_testVcfs.put(key, "org/pharmgkb/pharmcat/haplotype/IFNL3/rs12979860CC.vcf");
  }

  private PharmCAT m_pharmcat;
  private Path m_outputDir;

  public static void main(String[] args) {
    CliHelper cliHelper = new CliHelper(MethodHandles.lookup().lookupClass())
        .addOption("o", "output-dir", "directory to output to", true, "o")
        .addOption("n", "annotation-dir", "directory of guideline annotations (JSON files)", true, "n");

    try {
      if (!cliHelper.parse(args)) {
        System.exit(1);
      }

      Path outputDir = cliHelper.getValidDirectory("o", true);
      Path annoDir = cliHelper.getValidDirectory("n", false);

      PipelineTest piplelineTest = new PipelineTest(outputDir, annoDir);
      piplelineTest.execute();
    } catch (Exception e) {
      e.printStackTrace();
    }
  }


  private PipelineTest(Path outputDir, Path annoDir) throws IOException {
    m_outputDir = outputDir;
    m_pharmcat = new PharmCAT(outputDir, null, annoDir).keepMatcherOutput();
  }

  private void execute() throws Exception {
    System.out.println("Run time: " + new Date());

    for (String key : sf_testVcfs.keySet()) {
      runPipeline(key, sf_astrolabe);
      System.out.println("Generated "+key);
    }
  }

  private void runPipeline(String fileRoot, Path astrolabePath) throws Exception {
    Path sampleVcf = writeVcf(m_outputDir.resolve(fileRoot+".vcf"), sf_testVcfs.get(fileRoot));
    m_pharmcat.execute(sampleVcf, astrolabePath, null);
  }

  private Path writeVcf(Path outputVcf, Collection<String> filesToInclude) {
    Path headerFile = PathUtils.getPathToResource(sf_headerFile);

    try (FileWriter writer = new FileWriter(outputVcf.toFile())) {
      Files.lines(headerFile).filter(l -> l.startsWith("#")).forEach(l -> {
        try {
          writer.write(l);
          writer.write("\n");
        } catch (IOException e) {
          throw new RuntimeException(e);
        }
      });

      for (String filepath : filesToInclude) {
        Files.lines(PathUtils.getPathToResource(filepath)).filter(l -> !l.startsWith("#")).forEach(l -> {
          try {
            writer.write(l);
            writer.write("\n");
          } catch (IOException e) {
            throw new RuntimeException(e);
          }
        });
      }
    } catch (IOException e) {
      e.printStackTrace();
    }

    return outputVcf;
  }
}

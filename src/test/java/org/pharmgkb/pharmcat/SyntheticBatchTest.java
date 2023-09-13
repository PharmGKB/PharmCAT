package org.pharmgkb.pharmcat;

import java.io.FileWriter;
import java.io.IOException;
import java.lang.invoke.MethodHandles;
import java.nio.file.Files;
import java.nio.file.Path;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.List;
import java.util.stream.Collectors;
import com.google.common.collect.Lists;
import org.checkerframework.checker.nullness.qual.Nullable;
import org.junit.platform.launcher.Launcher;
import org.junit.platform.launcher.LauncherDiscoveryRequest;
import org.junit.platform.launcher.core.LauncherDiscoveryRequestBuilder;
import org.junit.platform.launcher.core.LauncherFactory;
import org.junit.platform.launcher.listeners.SummaryGeneratingListener;
import org.pharmgkb.common.util.CliHelper;
import org.pharmgkb.common.util.PathUtils;
import org.pharmgkb.pharmcat.reporter.model.DataSource;
import org.pharmgkb.pharmcat.util.CliUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import static org.junit.platform.engine.discovery.DiscoverySelectors.selectClass;


/**
 * This test generates a collection of synthetic sample VCF files from our test resources and then runs the matcher and
 * reporter on those samples, writing the output to a desired directory.
 * <p>
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
  private final boolean m_compact;
  private final List<DataSource> m_sources;


  public static void main(String[] args) {
    CliHelper cliHelper = new CliHelper(MethodHandles.lookup().lookupClass())
        .addOption("o", "output-dir", "directory to output to", false, "o")
        .addOption("a", "all-tests", "run all tests for Katrin")
        .addOption("re", "reporter-extended", "output extended report")
        .addOption("cpic", "cpic", "CPIC reports")
        .addOption("dpwg", "dpwg", "DPWG reports")
        .addOption("mega", "mega", "generate all variations in one run")
        ;

    try {
      if (!cliHelper.parse(args)) {
        System.exit(1);
      }

      if (cliHelper.hasOption("mega")) {
        System.out.println("Doing MEGA test set!");
        Path dir = null;
        if (cliHelper.hasOption("o")) {
          dir = cliHelper.getValidDirectory("o", true);
        } else {
          System.out.println("Must specify output directory (-o) when using -mega flag");
          System.exit(1);
        }

        List<DataSource> sources = Lists.newArrayList(DataSource.CPIC, DataSource.DPWG);
        doRun(dir.resolve("default"), true, sources, true);
        doRun(dir.resolve("extended"), false, sources, true);
        doRun(dir.resolve("cpic"), true, Lists.newArrayList(DataSource.CPIC), true);
        doRun(dir.resolve("cpic-extended"), false, Lists.newArrayList(DataSource.CPIC), true);
        doRun(dir.resolve("dpwg"), true, Lists.newArrayList(DataSource.DPWG), true);
        doRun(dir.resolve("dpwg-extended"), false, Lists.newArrayList(DataSource.DPWG), true);


      } else {
        List<DataSource> sources = Lists.newArrayList(DataSource.CPIC, DataSource.DPWG);
        boolean cpic = cliHelper.hasOption("cpic");
        boolean dpwg = cliHelper.hasOption("dpwg");
        if (cpic || dpwg) {
          sources.clear();
          if (cpic) {
            sources.add(DataSource.CPIC);
          }
          if (dpwg) {
            sources.add(DataSource.DPWG);
          }
        }
        Path dir = null;
        if (cliHelper.hasOption("o")) {
          dir = cliHelper.getValidDirectory("o", true);
        }
        doRun(dir, !cliHelper.hasOption("re"), sources, cliHelper.hasOption("a"));
      }

    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  private static void doRun(@Nullable Path dir, boolean compact, List<DataSource> sources, boolean allTests)
      throws Exception {
    if (dir != null) {
      System.out.println("Saving results to " + dir);
      TestUtils.setTestOutputDir(dir);
    }
    TestUtils.setSaveTestOutput(true);

    SyntheticBatchTest synthBatchTest = new SyntheticBatchTest(compact, sources);
    synthBatchTest.execute();

    if (allTests) {
      PipelineWrapper.setCompact(compact);
      PipelineWrapper.setSources(sources);
      SummaryGeneratingListener listener = new SummaryGeneratingListener();
      LauncherDiscoveryRequest request = LauncherDiscoveryRequestBuilder.request()
          .selectors(
              selectClass(PharmCATTest.class),
              selectClass(PipelineTest.class)
          )
          .build();
      Launcher launcher = LauncherFactory.create();
      launcher.discover(request);
      launcher.registerTestExecutionListeners(listener);
      launcher.execute(request);
    }
  }



  private void execute() throws Exception {

    makeReport("example.CYP2C19_noCYP2D6", new String[]{
        "cyp2c19/s1s1.vcf"
    }, null);

    makeReport("slco1b1.5.15", new String[]{
        "SLCO1B1/s5s15.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("slco1b1.s1s1", new String[]{
        "SLCO1B1/s1s1.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("slco1b1.s1s44", new String[]{
        "SLCO1B1/s1s44.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("slco1b1.s2s2", new String[]{
        "SLCO1B1/s2s2.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("slco1b1.s14_s19_t1", new String[]{
        "SLCO1B1/SLCO1B1_s14_s19_t1.vcf"
    }, null);

    makeReport("slco1b1.s15_s45_t1", new String[]{
        "SLCO1B1/SLCO1B1_s15_s45_t1.vcf"
    }, null);

    makeReport("slco1b1.s5s45", new String[]{
        "SLCO1B1/s5s45.vcf"
    }, null);

    makeReport("slco1b1.s19_s37_t1", new String[]{
        "SLCO1B1/SLCO1B1_s19_s37_t1.vcf"
    }, null);

    makeReport("slco1b1.s1_s14_t1", new String[]{
        "SLCO1B1/SLCO1B1_s1_s14_t1.vcf"
    }, null);

    makeReport("slco1b1.s1_s15_t1", new String[]{
        "SLCO1B1/SLCO1B1_s1_s15_t1.vcf"
    }, null);

    makeReport("slco1b1.s1_s25_t1", new String[]{
        "SLCO1B1/SLCO1B1_s1_s25_t1.vcf"
    }, null);

    makeReport("slco1b1.s37_s40_t1", new String[]{
        "SLCO1B1/SLCO1B1_s37_s40_t1.vcf"
    }, null);

    makeReport("slco1b1.s5_s14_t1", new String[]{
        "SLCO1B1/SLCO1B1_s5_s14_t1.vcf"
    }, null);

    makeReport("slco1b1.multi", new String[]{
        "SLCO1B1/multi.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("cyp2c19.onlyRs12769205", new String[]{
        "cyp2c19/rs12769205only.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("cyp2c19.refRs12769205", new String[]{
        "cyp2c19/rs12769205ref.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("cyp2c19.s2s3", new String[]{
        "cyp2c19/s2s3.vcf"
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
        "cyp2c19/s1s35.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("tpmt.star1s", new String[]{
        "TPMT/s1ss1ss3.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("tpmt.s1s1s", new String[]{
        "TPMT/s1s1s.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("tpmt.hom1s_het3a", new String[]{
        "TPMT/hom1s_het3a.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("tpmt.het3a", new String[]{
        "TPMT/het3a.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("tpmt.s15offdata", new String[]{
        "TPMT/s15offdata.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("tpmt.unexpected.allele", new String[]{
        "TPMT/unexpected.allele.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("tpmt.unexpected.allele.with.s9", new String[]{
        "TPMT/unexpected.allele.with.s9.vcf"
    }, sf_outsideCYP2D6File);

    makeReport("cyp2c19.rs28399504missing", new String[]{
        "cyp2c19/s4bs17rs28399504missing.vcf"
    }, sf_outsideCYP2D6File);


    makeReport("dpyd.stars12b", new String[]{
        "DPYD/s1s2b.vcf"
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

    makeReport("cyp2c9.s2s24", new String[]{
        "cyp2c9/s2s24.vcf"
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

    makeReportWithOutputString("g6pd.singleB", new String[]{
        "CYP2C9/s1s1.vcf"
    }, "G6PD\tB (wildtype)");

    makeReport("g6pd.homB", new String[]{
        "CYP2C9/s1s1.vcf"
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

    makeReportWithOutputString("cyp2d6.multiple", new String[]{
        "cyp2c9/s1s61.vcf"
    }, "CYP2D6\t*1/*1\nCYP2D6\t*1/*2\n");

    makeReportWithOutputString("mtrnr1.increased", new String[]{
        "cyp2c9/s1s61.vcf"
    }, "MT-RNR1\t1555A>G\n");

    sf_logger.info("Wrote reports to {}", TestUtils.getTestOutputDir());
  }


  private SyntheticBatchTest(boolean compact, List<DataSource> sources) throws IOException {
    m_compact = compact;
    m_sources = sources;

    String readmeContent = String.format("""
        # PharmCAT Example Reports
        
        Generated on: %s
        PharmCAT Version: %s",
        Sources: %s
        Style: %s
            """, new SimpleDateFormat("MMMMM dd, yyyy").format(new Date()), CliUtils.getVersion(),
        (compact ? "Compact" : "Full"),
        sources.stream().map(DataSource::toString).collect(Collectors.joining(", "))

    );
    Files.writeString(TestUtils.createTestFile(getClass(), "README.md"), readmeContent);
  }

  private void makeReport(String key, String[] testVcfs, Path outsideCallPath) throws Exception {
    Path testDir = TestUtils.getTestOutputDir(getClass(), false);
    if (!Files.exists(testDir)) {
      if (!testDir.toFile().mkdirs()) {
        throw new RuntimeException("Output directory could not be created " + testDir.toAbsolutePath());
      }
    }

    Path sampleVcf = writeVcf(testDir.resolve(key + ".vcf"), testVcfs);
    new Pipeline(new Env(),
        true, new VcfFile(sampleVcf), null, true,
        true, false, false, true,
        true, null, outsideCallPath,
        true, null, null, m_sources, m_compact, false, true,
        testDir, null, m_compact,
        Pipeline.Mode.TEST, null, false
    ).call();
  }

  private void makeReportWithOutputString(String key, String[] testVcfs, String outsideCalls) throws Exception {
    Path outsideCallPath = TestUtils.createTestFile(getClass(), "outsideCall.tsv");
    try (FileWriter fw = new FileWriter(outsideCallPath.toFile())) {
      fw.write(outsideCalls);
    }
    makeReport(key, testVcfs, outsideCallPath);
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

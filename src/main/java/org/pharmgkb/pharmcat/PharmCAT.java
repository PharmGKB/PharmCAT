package org.pharmgkb.pharmcat;

import java.lang.invoke.MethodHandles;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Objects;
import com.google.common.base.Stopwatch;
import org.pharmgkb.common.util.CliHelper;
import org.pharmgkb.common.util.TimeUtils;
import org.pharmgkb.pharmcat.util.CliUtils;


/**
 * Class to run the PharmCAT tool from input VCF file to final output report.
 *
 * @author Ryan Whaley
 */
public class PharmCAT {

  public static void main(String[] args) {
    Stopwatch stopwatch = Stopwatch.createStarted();

    try {
      CliHelper cliHelper = new CliHelper(MethodHandles.lookup().lookupClass())
          .addVersion("PharmCAT " + CliUtils.getVersion())
          // inputs
          .addOption("s", "samples", "Comma-separated list of samples", false, "samples")
          .addOption("S", "sample-file", "File containing a list of sample, one per line", false, "file")

          // named allele matcher args
          .addOption("matcher", "matcher", "Run named allele matcher independently")
          .addOption("vcf", "matcher-vcf", "Input VCF file for named allele matcher", false, "file")
          .addOption("ma", "matcher-all-results", "Return all possible diplotypes, not just top hits")
          .addOption("matcherHtml", "matcher-save-html", "Save named allele matcher results as HTML")

          // phenotyper args
          .addOption("phenotyper", "phenotyper", "Run phenotyper independently")
          .addOption("pi", "phenotyper-input", "JSON results from named allele matcher", false, "file")
          .addOption("po", "phenotyper-outside-call-file", "Path to an outside call file (TSV)", false, "file")

          // reporter args
          .addOption("reporter", "reporter", "Run reporter independently")
          .addOption("ri", "reporter-input", "JSON results from phenotyper", false, "file")
          .addOption("rt", "reporter-title", "Text to add to the report title", false, "title")
          .addOption("rs", "reporter-sources", "Comma-separated list of sources to limit report to: [CPIC, DPWG]", false, "sources")
          .addOption("re", "reporter-extended", "Output extended report")
          .addOption("reporterJson", "reporter-save-json", "Save reporter results as JSON")

          // outputs
          .addOption("o", "output-dir", "Directory to output to (optional, default is input file directory)", false, "directory")
          .addOption("bf", "base-filename", "The base name (without file extensions) used for output files, will default to base filename of input if not specified", false, "name")
          .addOption("del", "delete-intermediate-files", "Delete intermediate output files")
          // controls
          .addOption("def", "definitions-dir", "Directory containing named allele definitions (JSON files)", false, "dir")
          .addOption("research", "research-mode", "Comma-separated list of research features to enable: [cyp2d6, combinations]", false, "type");
      if (!cliHelper.parse(args)) {
        failIfNotTest();
        return;
      }

      BaseConfig config = new BaseConfig(cliHelper);

      VcfFile vcfFile = null;
      if (config.runMatcher) {
        if (cliHelper.hasOption("vcf")) {
          vcfFile = new VcfFile(cliHelper.getValidFile("vcf", true));
        } else {
          System.out.println(
              """
                  No input for Named Allele Matcher!

                  Please specify a VCF file (-vcf)"""
          );
          failIfNotTest();
          return;
        }
        if (cliHelper.hasOption("pi")) {
          System.out.println("Cannot specify phenotyper-input (-pi) if running named allele matcher");
          failIfNotTest();
          return;
        }
      } else if (!config.samples.isEmpty()) {
        throw new ReportableException("Cannot specify samples unless running matcher.");
      }

      Path phenotyperInputFile = null;
      Path phenotyperOutsideCallsFile = null;
      if (config.runPhenotyper) {
        if (cliHelper.hasOption("pi")) {
          phenotyperInputFile = cliHelper.getValidFile("pi", true);
        }
        if (cliHelper.hasOption("po")) {
          phenotyperOutsideCallsFile = cliHelper.getValidFile("po", true);
        }

        if (vcfFile == null && phenotyperInputFile == null && phenotyperOutsideCallsFile == null) {
          System.out.println("""
              No input for Phenotyper!

              Either:
                1. Run named allele matcher with VCF input, or
                2. Specify phenotyper-input (-pi) and/or phenotyper-outside-call-file (-po)"""
          );
          failIfNotTest();
          return;
        }
      }

      Path reporterInputFile = null;
      if (config.runReporter) {
        if (cliHelper.hasOption("ri")) {
          reporterInputFile = cliHelper.getValidFile("ri", true);
        }

        if (vcfFile == null && phenotyperInputFile == null && phenotyperOutsideCallsFile == null &&
            reporterInputFile == null) {
          System.out.println(
              """
                  No input for Reporter!

                  Either:
                    1. Run phenotyper, or
                    2. Specify reporter-input (-ri)"""
          );
          failIfNotTest();
          return;
        }
      }

      Env env = new Env(config.definitionDir);

      if (config.runMatcher) {
        Objects.requireNonNull(vcfFile);
        if (config.samples.isEmpty()) {
          config.samples.addAll(vcfFile.getSamples());
        } else {
          List<String> missing = new ArrayList<>();
          for (String sid : config.samples) {
            if (!vcfFile.getSamples().contains(sid)) {
              missing.add(sid);
            }
          }
          if (!missing.isEmpty()) {
            throw new ReportableException("The following samples could not be found in the VCF file: " +
                String.join(", ", missing));
          }
        }

        List<String> blankRuns = new ArrayList<>();
        int x = 0;
        boolean singleSample = config.samples.size() == 1;
        for (String sampleId : config.samples) {
          x += 1;
          if (config.samples.size() > 1) {
            System.out.println(x + " / " + config.samples.size() + " - " + sampleId);
          }
          Pipeline pipeline = new Pipeline(env,
              config.runMatcher, vcfFile, sampleId, singleSample,
              config.topCandidateOnly, config.callCyp2d6, config.findCombinations, config.matcherHtml,
              config.runPhenotyper, phenotyperInputFile, phenotyperOutsideCallsFile,
              config.runReporter, reporterInputFile, config.reporterTitle,
              config.reporterSources, config.reporterCompact, config.reporterJson, config.reporterHtml,
              config.outputDir, config.baseFilename, config.deleteIntermediateFiles,
              Pipeline.Mode.CLI, null, cliHelper.isVerbose());
          if (pipeline.call().getStatus() == PipelineResult.Status.NOOP) {
            failIfNotTest();
            blankRuns.add(sampleId);
          }

          if (x != config.samples.size() && config.samples.size() > 1) {
            System.out.println();
            System.out.println("---");
            System.out.println();
          }
        }
        if (!blankRuns.isEmpty()) {
          System.out.println("Nothing to do for " + String.join(", ", blankRuns));
        }

      } else {
        Pipeline pipeline = new Pipeline(env,
            false, null, null, true,
            config.topCandidateOnly, config.callCyp2d6, config.findCombinations, config.matcherHtml,
            config.runPhenotyper, phenotyperInputFile, phenotyperOutsideCallsFile,
            config.runReporter, reporterInputFile, config.reporterTitle,
            config.reporterSources, config.reporterCompact, config.reporterJson, config.reporterHtml,
            config.outputDir, config.baseFilename, config.deleteIntermediateFiles,
            Pipeline.Mode.CLI, null, cliHelper.isVerbose());
        if (pipeline.call().getStatus() == PipelineResult.Status.NOOP) {
          cliHelper.printHelp();
          System.out.println("Nothing to do.");
          failIfNotTest();
          return;
        }
      }

      System.out.println("Done.");
      if (cliHelper.isVerbose()) {
        System.out.println("Elapsed time: " + TimeUtils.humanReadablePreciseDuration(stopwatch.elapsed()));
      }

    } catch (CliHelper.InvalidPathException | ReportableException ex) {
      System.out.println(ex.getMessage());
      failIfNotTest();
    } catch (Exception e) {
      //noinspection CallToPrintStackTrace
      e.printStackTrace();
      failIfNotTest();
    }
  }


  /**
   * Only use {@link System#exit(int)} if not running from within test.
   */
  static void failIfNotTest() {
    try {
      Class.forName("org.pharmgkb.pharmcat.TestUtils");
    } catch (Exception ex) {
      System.exit(1);
    }
  }
}

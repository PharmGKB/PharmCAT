package org.pharmgkb.pharmcat;

import java.lang.invoke.MethodHandles;
import java.nio.file.Path;
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
          // named allele matcher args
          .addOption("matcher", "matcher", "run named allele matcher")
          .addOption("vcf", "matcher-vcf", "input VCF file for named allele matcher", false, "file")
          .addOption("ma", "matcher-all-results", "return all possible diplotypes, not just top hits")
          .addOption("matcherHtml", "matcher-save-html", "save named allele matcher results as HTML")

          // phenotyper args
          .addOption("phenotyper", "phenotyper", "run phenotyper")
          .addOption("pi", "phenotyper-input", "JSON results from named allele matcher", false, "file")
          .addOption("po", "phenotyper-outside-call-file", "path to an outside call file (TSV)", false, "file")

          // reporter args
          .addOption("reporter", "reporter", "run reporter")
          .addOption("ri", "reporter-input", "JSON results from phenotyper", false, "file")
          .addOption("rt", "reporter-title", "optional, text to add to the report title", false, "title")
          .addOption("rs", "reporter-sources", "comma-separated list of sources to limit recommendations to", false, "sources")
          .addOption("re", "reporter-extended", "generate extended report")
          .addOption("reporterJson", "reporter-save-json", "save reporter results as JSON")

          // outputs
          .addOption("o", "output-dir", "directory to output to (optional, default is input file directory)", false, "directory")
          .addOption("bf", "base-filename", "the base name (without file extensions) used for output files, will default to base filename of input if not specified", false, "name")
          // controls
          .addOption("def", "definitions-dir", "directory containing named allele definitions (JSON files)", false, "dir")
          .addOption("del", "delete-intermediary-files", "delete intermediary output files")
          .addOption("research", "research-mode", "comma-separated list of research features to enable [cyp2d6, combinations]", false, "type");
      if (!cliHelper.parse(args)) {
        failIfNotTest();
        return;
      }

      BaseConfig config = new BaseConfig(cliHelper);

      Path vcfFile = null;
      if (config.runMatcher) {
        if (cliHelper.hasOption("vcf")) {
          vcfFile = cliHelper.getValidFile("vcf", true);
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

      Pipeline pipeline = new Pipeline(new Env(config.definitionDir),
          config.runMatcher, vcfFile, config.topCandidateOnly, config.callCyp2d6, config.findCombinations,
          config.matcherHtml,
          config.runPhenotyper, phenotyperInputFile, phenotyperOutsideCallsFile,
          config.runReporter, reporterInputFile, config.reporterTitle,
          config.reporterSources, config.reporterCompact, config.reporterJson,
          config.outputDir, config.baseFilename, config.deleteIntermediateFiles, Pipeline.Mode.CLI);

      if (!pipeline.call()) {
        cliHelper.printHelp();
        failIfNotTest();
        return;
      }

      System.out.println("Done.");
      if (cliHelper.isVerbose()) {
        System.out.println("Elapsed time: " + TimeUtils.humanReadablePreciseDuration(stopwatch.elapsed()));
      }

    } catch (CliHelper.InvalidPathException | ReportableException ex) {
      System.out.println(ex.getMessage());
      failIfNotTest();
    } catch (Exception e) {
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

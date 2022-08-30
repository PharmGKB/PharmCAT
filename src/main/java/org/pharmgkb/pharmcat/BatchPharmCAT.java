package org.pharmgkb.pharmcat;

import java.lang.invoke.MethodHandles;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;
import java.util.TreeSet;
import java.util.stream.Stream;
import org.pharmgkb.common.util.CliHelper;
import org.pharmgkb.pharmcat.util.CliUtils;


/**
 * Command-line tool to run PharmCAT in batch mode.
 *
 * @author Mark Woon
 */
public class BatchPharmCAT {


  public static void main(String[] args) {
    try {
      CliHelper cliHelper = new CliHelper(MethodHandles.lookup().lookupClass())
          .addVersion("PharmCAT " + CliUtils.getVersion())
          // named allele matcher args
          .addOption("d", "dir", "directory containing source data files", true, "dir")

          .addOption("matcher", "matcher", "run named allele matcher")
          .addOption("ma", "matcher-all-results", "return all possible diplotypes, not just top hits")
          .addOption("md", "matcher-definitions-dir", "directory containing named allele definitions (JSON files)", false, "dir")
          .addOption("matcherHtml", "matcher-save-html", "save named allele matcher results as HTML")

          // phenotyper args
          .addOption("phenotyper", "phenotyper", "run phenotyper")

          // reporter args
          .addOption("reporter", "reporter", "run reporter")
          .addOption("rt", "reporter-title", "optional, text to add to the report title", false, "title")
          .addOption("rs", "reporter-sources", "comma-separated list of sources to limit report to", false, "sources")
          .addOption("rc", "reporter-compact", "output compact report")
          .addOption("reporterJson", "reporter-save-json", "save reporter results as JSON")

          // outputs
          .addOption("o", "output-dir", "directory to output to (optional, default is input file directory)", false, "directory")
          .addOption("bf", "base-filename", "the base name (without file extensions) used for output files, will default to base filename of input if not specified", false, "name")
          // controls
          .addOption("del", "delete-intermediary-files", "delete intermediary output files")
          .addOption("research", "research-mode", "comma-separated list of research features to enable [cyp2d6, combinations]", false, "type");
      if (!cliHelper.parse(args)) {
        System.exit(1);
      }

      Path dir = cliHelper.getValidDirectory("dir", false);
      PharmCAT.BaseConfig config = new PharmCAT.BaseConfig(cliHelper);

      Collection<String> keys = null;

      Map<String, Path> vcfFiles = new HashMap<>();
      if (config.runMatcher) {
        try (Stream<Path> stream = Files.list(dir)) {
          stream.filter(f -> f.toString().toLowerCase().endsWith(".vcf"))
              .forEach(f -> {
                String baseFilename = PharmCAT.getBaseFilename(f);
                vcfFiles.put(baseFilename, f);
              });
        }
        if (vcfFiles.size() == 0) {
          System.out.println("No VCF files found in " + dir);
          return;
        }
        keys = vcfFiles.keySet();
        System.out.println("Found " + keys.size() + " VCF files...");
      }

      Map<String, Path> phenotyperInputFiles = new HashMap<>();
      Map<String, Path> phenotyperOutsideCallFiles = new HashMap<>();
      if (config.runPhenotyper) {
        if (keys == null) {
          try (Stream<Path> stream = Files.list(dir)) {
            stream.filter(f -> f.toString().toLowerCase().endsWith(".match.json"))
                .forEach(f -> {
                  String baseFilename = PharmCAT.getBaseFilename(f);
                  phenotyperInputFiles.put(baseFilename, f);
                });
          }
          if (phenotyperInputFiles.size() == 0) {
            System.out.println("No Phenotyper input files found in " + dir);
            return;
          }
          keys = phenotyperInputFiles.keySet();
          System.out.println("Found " + keys.size() + " Phenotyper input files...");
        }

        for (String key : keys) {
          Path f = dir.resolve(key + ".outside.tsv");
          if (Files.isRegularFile(f)) {
            phenotyperOutsideCallFiles.put(key, f);
          }
        }
        if (phenotyperOutsideCallFiles.size() > 0) {
          System.out.println("Found " + keys.size() + " outside call files...");
        }
      }

      Map<String, Path> reporterInputFiles = new HashMap<>();
      if (config.runReporter) {
        if (keys == null) {
          try (Stream<Path> stream = Files.list(dir)) {
            stream.filter(f -> f.toString().toLowerCase().endsWith(".phenotype.json"))
                .forEach(f -> {
                  String baseFilename = PharmCAT.getBaseFilename(f);
                  reporterInputFiles.put(baseFilename, f);
                });
          }
          if (reporterInputFiles.size() == 0) {
            System.out.println("No Reporter input files found in " + dir);
            return;
          }
          keys = reporterInputFiles.keySet();
          System.out.println("Found " + keys.size() + " Reporter input files...");
        }
      }

      if (keys == null || keys.size() == 0) {
        System.out.println("No files found to process.");
        return;
      }

      int x = 0;
      for (String key : new TreeSet<>(keys)) {
        Path vcfFile = vcfFiles.get(key);
        Path phenotyperInputFile = phenotyperInputFiles.get(key);
        Path phenotyperOutsideCallsFile = phenotyperOutsideCallFiles.get(key);
        Path reporterInputFile = reporterInputFiles.get(key);

        x += 1;
        System.out.println(x + ": " + key);
        new PharmCAT(config.runMatcher, vcfFile, config.definitionReader,
            config.topCandidateOnly, config.callCyp2d6, config.findCombinations, config.matcherHtml,
            config.runPhenotyper, phenotyperInputFile, phenotyperOutsideCallsFile,
            config.runReporter, reporterInputFile, config.reporterTitle,
            config.reporterSources, config.reporterCompact, config.reporterJson,
            config.outputDir, config.baseFilename, config.deleteIntermediateFiles, PharmCAT.Mode.CLI)
            .execute();
        System.out.println("---");
      }

    } catch (CliHelper.InvalidPathException ex) {
      System.out.println(ex.getMessage());
      System.exit(1);
    } catch (Exception e) {
      e.printStackTrace();
      System.exit(1);
    }
  }
}

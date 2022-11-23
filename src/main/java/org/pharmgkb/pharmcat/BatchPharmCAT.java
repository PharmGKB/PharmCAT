package org.pharmgkb.pharmcat;

import java.lang.invoke.MethodHandles;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.stream.Stream;
import com.google.common.base.Stopwatch;
import org.pharmgkb.common.util.CliHelper;
import org.pharmgkb.common.util.TimeUtils;
import org.pharmgkb.pharmcat.haplotype.VcfSampleReader;
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
          .addOption("i", "input-dir", "directory containing source data files", true, "dir")

          .addOption("matcher", "matcher", "run named allele matcher")
          .addOption("ma", "matcher-all-results", "return all possible diplotypes, not just top hits")
          .addOption("md", "matcher-definitions-dir", "directory containing named allele definitions (JSON files)", false, "dir")
          .addOption("matcherHtml", "matcher-save-html", "save named allele matcher results as HTML")

          // phenotyper args
          .addOption("phenotyper", "phenotyper", "run phenotyper")

          // reporter args
          .addOption("reporter", "reporter", "run reporter")
          .addOption("rs", "reporter-sources", "comma-separated list of sources to limit report to", false, "sources")
          .addOption("re", "reporter-extended", "output extended report")
          .addOption("reporterJson", "reporter-save-json", "save reporter results as JSON")

          // outputs
          .addOption("o", "output-dir", "directory to output to (optional, default is input file directory)", false, "directory")
          .addOption("bf", "base-filename", "the base name (without file extensions) used for output files, will default to base filename of input if not specified", false, "name")
          // controls
          .addOption("cp", "max-concurrent-processes", "maximum number of processes to use", false, "num")
          .addOption("def", "definitions-dir", "directory containing named allele definitions (JSON files)", false, "dir")
          .addOption("del", "delete-intermediary-files", "delete intermediary output files")
          .addOption("research", "research-mode", "comma-separated list of research features to enable [cyp2d6, combinations]", false, "type");
      if (!cliHelper.parse(args)) {
        PharmCAT.failIfNotTest();
        return;
      }

      BaseConfig config = new BaseConfig(cliHelper);

      Path dir = cliHelper.getValidDirectory("i", false);
      List<Path> allVcfFiles = new ArrayList<>();
      List<Path> allMatchFiles = new ArrayList<>();
      List<Path> allOutsideCallFiles = new ArrayList<>();
      List<Path> allPhenotypeFiles = new ArrayList<>();
      try (Stream<Path> stream = Files.list(dir)) {
        stream.filter(Files::isRegularFile)
            .forEach(f -> {
              String name = f.toString().toLowerCase();
              if (name.endsWith(".vcf")) {
                allVcfFiles.add(f);
              } else if (name.endsWith(".match.json")) {
                allMatchFiles.add(f);
              } else if (name.endsWith(".outside.tsv")) {
                allOutsideCallFiles.add(f);
              } else if (name.endsWith(".phenotype.json")) {
                allPhenotypeFiles.add(f);
              }
            });
      }

      final SortedSet<String> keys = new TreeSet<>();
      Map<String, Path> vcfFiles = new HashMap<>();
      if (config.runMatcher) {
        allVcfFiles.forEach(f -> {
          String baseFilename = BaseConfig.getBaseFilename(f);
          vcfFiles.put(baseFilename, f);
        });
        if (vcfFiles.size() == 0) {
          System.out.println("No VCF files found in " + dir);
          return;
        }
        keys.addAll(vcfFiles.keySet());
        System.out.println("Found " + vcfFiles.size() + " VCF files...");
      }

      Map<String, Path> phenotyperInputFiles = new HashMap<>();
      Map<String, Path> phenotyperOutsideCallFiles = new HashMap<>();
      if (config.runPhenotyper) {
        if (keys.isEmpty()) {
          allMatchFiles.forEach(f -> {
            String baseFilename = BaseConfig.getBaseFilename(f);
            phenotyperInputFiles.put(baseFilename, f);
          });
          allOutsideCallFiles.forEach(f -> {
            String baseFilename = BaseConfig.getBaseFilename(f);
            if (!keys.contains(baseFilename)) {
              System.out.println("Warning: outside call file with no matching .match.json - " + f.getFileName());
            }
            phenotyperOutsideCallFiles.put(baseFilename, f);
          });
          if (phenotyperInputFiles.size() == 0 && phenotyperOutsideCallFiles.size() == 0) {
            System.out.println("No phenotyper input files/outside call files found in " + dir);
            return;
          }
          keys.addAll(phenotyperInputFiles.keySet());
          keys.addAll(phenotyperOutsideCallFiles.keySet());

        } else {
          allMatchFiles.forEach(f -> {
            String baseFilename = BaseConfig.getBaseFilename(f);
            if (keys.contains(baseFilename)) {
              System.out.println("Ignoring " + f.getFileName() + " - will recompute");
              return;
            }
            // allow mixing phenotyper-only inputs
            phenotyperInputFiles.put(baseFilename, f);
            keys.add(baseFilename);
          });
          allOutsideCallFiles.forEach(f -> {
            String baseFilename = BaseConfig.getBaseFilename(f);
            if (!keys.contains(baseFilename)) {
              System.out.println("Warning: outside call file with no matching .vcf or .match.json - " + f.getFileName());
              keys.add(baseFilename);
            }
            phenotyperOutsideCallFiles.put(baseFilename, f);
          });
        }

        if (phenotyperInputFiles.size() > 0) {
          System.out.println(getCountString(phenotyperInputFiles, allMatchFiles) + " phenotyper input files...");
        }
        if (phenotyperOutsideCallFiles.size() > 0) {
          System.out.println(getCountString(phenotyperOutsideCallFiles, allOutsideCallFiles) + " outside call files...");
        }
      }

      Map<String, Path> reporterInputFiles = new HashMap<>();
      if (config.runReporter) {
        if (keys.isEmpty()) {
          // reporter only
          allPhenotypeFiles.forEach(f -> {
            String baseFilename = BaseConfig.getBaseFilename(f);
            reporterInputFiles.put(baseFilename, f);
          });
          if (reporterInputFiles.size() == 0) {
            System.out.println("No reporter input files found in " + dir);
            return;
          }
          keys.addAll(reporterInputFiles.keySet());
          System.out.println(getCountString(reporterInputFiles, allPhenotypeFiles) + " reporter input files...");

        } else {
          allPhenotypeFiles.forEach(f -> {
            String baseFilename = BaseConfig.getBaseFilename(f);
            if (keys.contains(baseFilename)) {
              System.out.println("Ignoring " + f.getFileName() + " - will recompute");
              return;
            }
            // allow mixing reporter-only inputs
            System.out.println("Warning: reporter input file with no matching .vcf, .match.json or .outside.tsv - " +
                f.getFileName());
            reporterInputFiles.put(baseFilename, f);
            keys.add(baseFilename);
          });
        }
      }

      if (keys.isEmpty()) {
        System.out.println("No files found to process.");
        return;
      }
      SortedSet<String> sortedKeys = new TreeSet<>(keys);
      System.out.println();
      System.out.println("Queueing up " + sortedKeys.size() + " samples to process...");

      Env env = new Env(config.definitionDir);
      List<Pipeline> tasks = new ArrayList<>();
      for (String key : new TreeSet<>(keys)) {
        Path vcfFile = vcfFiles.get(key);
        Path phenotyperInputFile = phenotyperInputFiles.get(key);
        Path phenotyperOutsideCallsFile = phenotyperOutsideCallFiles.get(key);
        Path reporterInputFile = reporterInputFiles.get(key);

        if (cliHelper.isVerbose()) {
          System.out.println("  * " + key);
        }
        boolean runMatcher = config.runMatcher;
        if (runMatcher && vcfFile == null) {
          // allow mixing phenotyper-only inputs
          runMatcher = false;
        }
        boolean runPhenotyper = config.runPhenotyper;
        if (runPhenotyper && vcfFile == null && phenotyperInputFile == null && phenotyperOutsideCallsFile == null) {
          runPhenotyper = false;
        }

        if (runMatcher) {
          VcfSampleReader sampleReader = new VcfSampleReader(vcfFile);
          if (sampleReader.getSamples().size() > 1) {
            for (String sampleId : sampleReader.getSamples()) {
              tasks.add(new Pipeline(env,
                  true, vcfFile, sampleId,
                  config.topCandidateOnly, config.callCyp2d6, config.findCombinations, config.matcherHtml,
                  runPhenotyper, phenotyperInputFile, phenotyperOutsideCallsFile,
                  config.runReporter, reporterInputFile, config.reporterTitle,
                  config.reporterSources, config.reporterCompact, config.reporterJson,
                  config.outputDir, config.baseFilename, config.deleteIntermediateFiles, Pipeline.Mode.BATCH));
            }
          } else {
            tasks.add(new Pipeline(env,
                true, vcfFile, null,
                config.topCandidateOnly, config.callCyp2d6, config.findCombinations, config.matcherHtml,
                runPhenotyper, phenotyperInputFile, phenotyperOutsideCallsFile,
                config.runReporter, reporterInputFile, config.reporterTitle,
                config.reporterSources, config.reporterCompact, config.reporterJson,
                config.outputDir, config.baseFilename, config.deleteIntermediateFiles, Pipeline.Mode.BATCH));
          }

        } else {
          tasks.add(new Pipeline(env,
              false, null, null,
              config.topCandidateOnly, config.callCyp2d6, config.findCombinations, config.matcherHtml,
              runPhenotyper, phenotyperInputFile, phenotyperOutsideCallsFile,
              config.runReporter, reporterInputFile, config.reporterTitle,
              config.reporterSources, config.reporterCompact, config.reporterJson,
              config.outputDir, config.baseFilename, config.deleteIntermediateFiles, Pipeline.Mode.BATCH));
        }
      }

      int maxProcesses = Runtime.getRuntime().availableProcessors() - 2;
      if (cliHelper.hasOption("cp")) {
        try {
          maxProcesses = cliHelper.getIntValue("cp");
          if (maxProcesses > Runtime.getRuntime().availableProcessors()) {
            maxProcesses = Runtime.getRuntime().availableProcessors();
          }
        } catch (NumberFormatException ex) {
          System.out.println("\"" + cliHelper.getValue("cp") + "\" is not an integer.");
          return;
        }
      }
      if (maxProcesses < 1) {
        maxProcesses = 1;
      }
      System.out.println();
      System.out.println("Running PharmCAT in batch mode with a maximum of " + maxProcesses + " processes.");
      Stopwatch stopwatch = Stopwatch.createStarted();
      ExecutorService executor = Executors.newFixedThreadPool(maxProcesses);
      List<Future<Boolean>> futures = executor.invokeAll(tasks);
      executor.shutdown();

      // must iterate through in case of errors
      for (Future<Boolean> future : futures) {
        future.get();
      }

      System.out.println();
      System.out.println("Done.");
      System.out.println("Elapsed time: " + TimeUtils.humanReadablePreciseDuration(stopwatch.elapsed()));

    } catch (CliHelper.InvalidPathException ex) {
      System.out.println(ex.getMessage());
      System.exit(1);
    } catch (Exception e) {
      e.printStackTrace();
      System.exit(1);
    }
  }


  private static String getCountString(Map<String, Path> using, Collection<Path> found) {
    if (using.size() != found.size()) {
      return "Adding " + using.size() + "/" + found.size();
    }
    return "Found " + using.size();
  }
}

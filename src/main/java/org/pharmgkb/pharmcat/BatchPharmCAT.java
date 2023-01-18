package org.pharmgkb.pharmcat;

import java.io.File;
import java.io.IOException;
import java.lang.invoke.MethodHandles;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.TreeMap;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import com.google.common.base.Preconditions;
import com.google.common.base.Stopwatch;
import org.apache.commons.io.FileUtils;
import org.checkerframework.checker.nullness.qual.Nullable;
import org.pharmgkb.common.util.CliHelper;
import org.pharmgkb.common.util.TimeUtils;
import org.pharmgkb.pharmcat.util.CliUtils;


/**
 * Command-line tool to run PharmCAT in batch mode.
 *
 * @author Mark Woon
 */
public class BatchPharmCAT {
  private static final int sf_procsPerGb = 16;
  private static final long sf_bytesPerProcess = (1024 / sf_procsPerGb) * 1024 * 1024;
  private final BaseConfig m_config;
  private final boolean m_verbose;
  private final Map<String, VcfFile> m_vcfFilesToProcess = new TreeMap<>();
  private final Map<String, Path> m_matchFilesToProcess = new TreeMap<>();
  private final Map<String, Path> m_outsideCallFilesToProcess = new TreeMap<>();
  private final Map<String, Path> m_phenotypeFilesToProcess = new TreeMap<>();


  public static void main(String[] args) {
    try {
      CliHelper cliHelper = new CliHelper(MethodHandles.lookup().lookupClass())
          .addVersion("PharmCAT " + CliUtils.getVersion())
          // inputs
          .addOption("i", "input-dir", "directory containing source data files", false, "dir")
          .addOption("s", "samples", "comma-separated list of samples", false, "samples")

          // named allele matcher args
          .addOption("matcher", "matcher", "run named allele matcher independently")
          .addOption("vcf", "matcher-vcf", "input VCF file for named allele matcher", false, "file")
          .addOption("ma", "matcher-all-results", "return all possible diplotypes, not just top hits")
          .addOption("md", "matcher-definitions-dir", "directory containing named allele definitions (JSON files)", false, "dir")
          .addOption("matcherHtml", "matcher-save-html", "save named allele matcher results as HTML")

          // phenotyper args
          .addOption("phenotyper", "phenotyper", "run phenotyper independently")

          // reporter args
          .addOption("reporter", "reporter", "run reporter independently")
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
      int maxProcesses = Runtime.getRuntime().availableProcessors() - 2;
      if (cliHelper.hasOption("cp")) {
        try {
          int cp = cliHelper.getIntValue("cp");
          if (cp > Runtime.getRuntime().availableProcessors()) {
            maxProcesses = Runtime.getRuntime().availableProcessors();
            System.out.println("Warning: This system only has " + maxProcesses + " processors.");
            System.out.println("Limiting concurrent processes to " + maxProcesses + ".");
            int futureMax = Math.min(1, maxProcesses - 2);
            System.out.println("Recommending maximum of '-cp " + futureMax + "' in the future.");
          } else {
            maxProcesses = cp;
          }
        } catch (NumberFormatException ex) {
          System.out.println("\"" + cliHelper.getValue("cp") + "\" is not an integer.");
          PharmCAT.failIfNotTest();
          return;
        }
      }
      if (maxProcesses < 1) {
        maxProcesses = 1;
      }

      Path inputDir = null;
      if (cliHelper.hasOption("i")) {
        inputDir = cliHelper.getValidDirectory("i", false);
      }
      Path vcfFile = null;
      if (cliHelper.hasOption("vcf")) {
        vcfFile = cliHelper.getValidFile("vcf", true);
        if (inputDir == null) {
          inputDir = vcfFile.getParent();
        }
      }
      if (inputDir == null) {
        System.err.println("Missing input (specify with -i and/or -vcf");
        PharmCAT.failIfNotTest();
        return;
      }

      BatchPharmCAT pcat = new BatchPharmCAT(config, inputDir, vcfFile, cliHelper.isVerbose());
      pcat.execute(maxProcesses);

    } catch (CliHelper.InvalidPathException | ReportableException ex) {
      System.out.println(ex.getMessage());
      PharmCAT.failIfNotTest();
    } catch (Exception e) {
      e.printStackTrace();
      PharmCAT.failIfNotTest();
    }
  }


  private BatchPharmCAT(BaseConfig config, Path inputDir, @Nullable Path vcfFile, boolean verbose)
      throws IOException, ReportableException {
    m_config = config;
    m_verbose = verbose;

    for (File f : Objects.requireNonNull(inputDir.toFile().listFiles())) {
      Path file = f.toPath();
      String name = file.toString().toLowerCase();
      String basename = BaseConfig.getBaseFilename(file);
      if (VcfFile.isVcfFile(file)) {
        if (config.runMatcher) {
          m_vcfFilesToProcess.put(basename, new VcfFile(file));
        }
      } else if (name.endsWith(".match.json")) {
        if (config.runPhenotyper) {
          m_matchFilesToProcess.put(basename, file);
        }
      } else if (name.endsWith(".outside.tsv")) {
        if (config.runPhenotyper) {
          m_outsideCallFilesToProcess.put(basename, file);
        }
      } else if (name.endsWith(".phenotype.json")) {
        if (config.runReporter) {
          m_phenotypeFilesToProcess.put(basename, file);
        }
      }
    }
    if (vcfFile != null) {
      if (!Files.isRegularFile(vcfFile)) {
        System.err.println("Not a file: " + vcfFile);
        PharmCAT.failIfNotTest();
        return;
      }
      // input VCF file trumps other VCF files in inputDir
      m_vcfFilesToProcess.clear();
      m_vcfFilesToProcess.put(BaseConfig.getBaseFilename(vcfFile), new VcfFile(vcfFile));
    }

    if (m_vcfFilesToProcess.isEmpty() && m_matchFilesToProcess.isEmpty() && m_outsideCallFilesToProcess.isEmpty() &&
        m_phenotypeFilesToProcess.isEmpty()) {
      List<String> types = new ArrayList<>();
      if (config.runMatcher) {
        types.add("*.vcf");
      }
      if (config.runPhenotyper) {
        types.add("*.match.json");
        types.add("*.outside.tsv");
      }
      if (config.runReporter) {
        types.add("*.phenotype.json");
      }
      throw new ReportableException("No input files (" + String.join(", ", types) + ") found in " + inputDir);
    }
  }


  private void execute(int maxProcesses) throws ExecutionException, InterruptedException, IOException,
      ReportableException {

    List<Builder> taskBuilders = new ArrayList<>();

    if (m_config.runMatcher) {
      System.out.println("Found " + m_vcfFilesToProcess.size() + " VCF files...");
      for (String baseFilename : m_vcfFilesToProcess.keySet()) {
        VcfFile vcfFile = m_vcfFilesToProcess.get(baseFilename);
        if (vcfFile != null) {
          boolean singleSample = vcfFile.getSamples().size() == 1;
          for (String sampleId : vcfFile.getSamples()) {
            if (m_config.runSample(sampleId)) {
              taskBuilders.add(new Builder().fromMatcher(baseFilename, vcfFile, sampleId, singleSample));
            }
          }
        }
      }
    }

    if (m_config.runPhenotyper) {
      if (m_matchFilesToProcess.size() > 0) {
        System.out.println("Found " + m_matchFilesToProcess.size() + " independent phenotyper input file" +
            (m_matchFilesToProcess.size() > 1 ? "s" : "") + "...");
        for (String baseFilename : m_matchFilesToProcess.keySet()) {
          if (m_config.runSample(baseFilename)) {
            taskBuilders.add(new Builder().fromPhenotyper(baseFilename));
          }
        }
      }
      if (m_outsideCallFilesToProcess.size() > 0) {
        System.out.println("Found " + m_outsideCallFilesToProcess.size() + " independent phenotyper outside call file" +
            (m_outsideCallFilesToProcess.size() > 1 ? "s" : "") + "...");
        for (String baseFilename : m_outsideCallFilesToProcess.keySet()) {
          if (m_config.runSample(baseFilename)) {
            taskBuilders.add(new Builder().fromPhenotyper(baseFilename));
          }
        }
      }
    }

    if (m_config.runReporter) {
      if (m_phenotypeFilesToProcess.size() > 0) {
        System.out.println("Found " + m_phenotypeFilesToProcess.size() + " independent reporter input file" +
            (m_phenotypeFilesToProcess.size() > 1 ? "s" : "") + "...");
        for (String baseFilename : m_phenotypeFilesToProcess.keySet()) {
          if (m_config.runSample(baseFilename)) {
            taskBuilders.add(new Builder().fromReporter(baseFilename));
          }
        }
      }
    }

    System.out.println();
    System.out.println("Queueing up " + taskBuilders.size() + " samples to process...");
    Env env = new Env(m_config.definitionDir);
    List<Pipeline> tasks = new ArrayList<>();
    Map<String, Pipeline> taskMap = new HashMap<>();
    for (Builder builder : taskBuilders) {
      Pipeline pipeline = builder.build(env);
      tasks.add(pipeline);
    }

    int processes = Math.min(tasks.size(), maxProcesses);
    System.out.println();
    System.out.println("Running PharmCAT in batch mode with a maximum of " + processes + " processes.");
    if (processes > 1) {
      long maxMem = Runtime.getRuntime().maxMemory();
      long memPerProcess = maxMem / processes;
      if (memPerProcess < sf_bytesPerProcess) {
        System.out.println("Warning: PharmCAT only has access to " + FileUtils.byteCountToDisplaySize(maxMem));
        int gb = processes / sf_procsPerGb;
        String recMem;
        if (gb > 0 && processes % sf_procsPerGb == 0) {
          recMem = gb + " G";
        } else {
          recMem = ((gb * 1024) + ((processes % sf_procsPerGb) * 1024 / sf_procsPerGb)) + " M";
        }
        System.out.println("Recommend boosting memory to PharmCAT to " + recMem + "B (using -Xmx" +
            recMem.replace(" ", "") + ")");
        long recCp = Math.min(1, maxMem / sf_bytesPerProcess);
        System.out.println("Or running with " + recCp + " process" + (recCp == 1L ? "" : "es") +
            " max (using -cp " + recCp + ")");
      }
    }

    if (m_verbose) {
      System.out.println("Current memory usage: " + FileUtils.byteCountToDisplaySize(
          Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory()));
    }
    Stopwatch stopwatch = Stopwatch.createStarted();
    ExecutorService executor = Executors.newWorkStealingPool(processes);
    List<Future<PipelineResult>> futures = executor.invokeAll(tasks);
    executor.shutdown();

    // must iterate through in case of errors
    for (Future<PipelineResult> future : futures) {
      PipelineResult rez = future.get();
      if (rez.getStatus() == PipelineResult.Status.FAILURE) {
        System.out.println("FAILED " + (rez.getSampleId() == null ? rez.getBasename() : rez.getSampleId()));
      }
    }

    System.out.println();
    System.out.println("Done.");
    System.out.println("Elapsed time: " + TimeUtils.humanReadablePreciseDuration(stopwatch.elapsed()));

  }


  public class Builder {
    private String m_baseFilename;
    private boolean m_runMatcher;
    private VcfFile m_vcfFile;
    private String m_sampleId;
    private boolean m_runPhenotyper;
    private Path m_piFile;
    private Path m_poFile;
    private boolean m_runReporter;
    private Path m_riFile;
    private boolean m_singleSample;


    public Builder fromMatcher(String baseFilename, VcfFile file, @Nullable String sampleId, boolean singleSample) {
      Preconditions.checkState(m_config.runMatcher);
      Preconditions.checkNotNull(file);
      m_baseFilename = baseFilename;
      m_vcfFile = file;
      m_sampleId = sampleId;
      m_runMatcher = true;
      if (sampleId == null) {
        findPhenotyperFiles(baseFilename);
        findReporterFiles(baseFilename);
      } else {
        findPhenotyperFiles(sampleId);
        findReporterFiles(sampleId);
        findPhenotyperFiles(baseFilename + "." + sampleId);
        findReporterFiles(baseFilename + "." + sampleId);
      }
      m_singleSample = singleSample;
      return this;
    }


    public Builder fromPhenotyper(String baseFilename) {
      Preconditions.checkState(m_config.runPhenotyper);
      m_baseFilename = baseFilename;
      findPhenotyperFiles(baseFilename);
      findReporterFiles(baseFilename);
      m_singleSample = true;
      return this;
    }


    public Builder fromReporter(String baseFilename) {
      Preconditions.checkState(m_config.runReporter);
      m_baseFilename = baseFilename;
      findReporterFiles(baseFilename);
      m_singleSample = true;
      return this;
    }


    public Pipeline build(Env env) throws ReportableException {
      Preconditions.checkState(m_runMatcher || m_runPhenotyper || m_runReporter);

      return new Pipeline(env,
          m_runMatcher, m_vcfFile, m_sampleId,
          m_config.topCandidateOnly, m_config.callCyp2d6, m_config.findCombinations, m_config.matcherHtml,
          m_runPhenotyper, m_piFile, m_poFile,
          m_runReporter, m_riFile, m_config.reporterTitle,
          m_config.reporterSources, m_config.reporterCompact, m_config.reporterJson,
          m_config.outputDir, m_config.baseFilename, m_config.deleteIntermediateFiles,
          Pipeline.Mode.BATCH, m_singleSample, m_verbose);
    }


    private void findPhenotyperFiles(String basename) {
      if (!m_config.runPhenotyper) {
        return;
      }
      // pi file
      if (m_matchFilesToProcess.containsKey(basename)) {
        Path file = m_matchFilesToProcess.get(basename);
        if (m_config.runMatcher && m_vcfFile != null) {
          System.out.println("Ignoring " + file.getFileName() + " - will recompute");
        } else {
          m_piFile = pickFirstFile(m_piFile, file);
        }
        m_matchFilesToProcess.remove(basename);
      }
      // po file
      if (m_outsideCallFilesToProcess.containsKey(basename)) {
        Path file = m_outsideCallFilesToProcess.get(basename);
        m_poFile = pickFirstFile(m_poFile, file);
        m_outsideCallFilesToProcess.remove(basename);

        if ((m_piFile == null && !m_config.runMatcher) || (m_piFile == null && m_vcfFile == null)) {
          System.out.println("Warning: lone outside call file (" + m_poFile.getFileName() +
              ") with no matching .vcf or .match.json");
        }
      }
      m_runPhenotyper = m_vcfFile != null || m_piFile != null || m_poFile != null;
    }

    private void findReporterFiles(String basename) {
      if (!m_config.runReporter) {
        return;
      }
      // ri file
      if (m_phenotypeFilesToProcess.containsKey(basename)) {
        Path file = m_phenotypeFilesToProcess.get(basename);
        if ((m_config.runMatcher && m_vcfFile != null) ||
            (m_config.runPhenotyper && (m_piFile != null || m_poFile != null))) {
          System.out.println("Ignoring " + file.getFileName() + " - will recompute");
        } else {
          m_riFile = pickFirstFile(m_riFile, file);
        }
        m_phenotypeFilesToProcess.remove(basename);
      }
      m_runReporter = m_vcfFile != null || m_piFile != null || m_poFile != null || m_riFile != null;
    }

    private Path pickFirstFile(Path origFile, Path newFile) {
      if (origFile == null) {
        return newFile;
      } else {
        System.out.println("Ignoring " + newFile.getFileName() + " - using " + origFile.getFileName() + " instead");
        return origFile;
      }
    }
  }
}

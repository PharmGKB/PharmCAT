package org.pharmgkb.pharmcat;

import java.io.IOException;
import java.lang.invoke.MethodHandles;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.concurrent.TimeUnit;
import java.util.stream.Collectors;
import com.google.common.base.Splitter;
import com.google.common.base.Stopwatch;
import org.apache.commons.io.FilenameUtils;
import org.checkerframework.checker.nullness.qual.Nullable;
import org.pharmgkb.common.util.CliHelper;
import org.pharmgkb.pharmcat.haplotype.DefinitionReader;
import org.pharmgkb.pharmcat.haplotype.NamedAlleleMatcher;
import org.pharmgkb.pharmcat.haplotype.ResultSerializer;
import org.pharmgkb.pharmcat.haplotype.model.GeneCall;
import org.pharmgkb.pharmcat.haplotype.model.Result;
import org.pharmgkb.pharmcat.phenotype.Phenotyper;
import org.pharmgkb.pharmcat.reporter.ReportContext;
import org.pharmgkb.pharmcat.reporter.format.HtmlFormat;
import org.pharmgkb.pharmcat.reporter.format.JsonFormat;
import org.pharmgkb.pharmcat.reporter.io.OutsideCallParser;
import org.pharmgkb.pharmcat.reporter.model.DataSource;
import org.pharmgkb.pharmcat.reporter.model.MessageAnnotation;
import org.pharmgkb.pharmcat.reporter.model.OutsideCall;
import org.pharmgkb.pharmcat.util.CliUtils;
import org.pharmgkb.pharmcat.util.DataManager;


/**
 * Class to run the PharmCAT tool from input VCF file to final output report.
 *
 * @author Ryan Whaley
 */
public class PharmCAT {
  private static final Splitter sf_commaSplitter = Splitter.on(",").trimResults().omitEmptyStrings();
  public enum Mode {
    /**
     * Default mode.  Prints informative messages to console.
     */
    CLI,
    /**
     * In test mode, PharmCAT tries not to include version/timestamp in output to simplify diffing.
     */
    TEST
  };
  private DefinitionReader m_definitionReader;
  private final boolean m_runMatcher;
  private Path m_vcfFile;
  private boolean m_topCandidateOnly = true;
  private boolean m_findCombinations;
  private boolean m_callCyp2d6;
  private Path m_matcherJsonFile;
  private Path m_matcherHtmlFile;

  private final boolean m_runPhenotyper;
  private Path m_phenotyperInputFile;
  private Path m_phenotyperOutsideCallsFile;
  private Path m_phenotyperJsonFile;

  private final boolean m_runReporter;
  private Path m_reporterInputFile;
  private String m_reporterTitle;
  private boolean m_reporterCompact;
  private List<DataSource> m_reporterSources;
  private Path m_reporterJsonFile;
  private Path m_reporterHtmlFile;
  private ReportContext m_reportContext;

  private final boolean m_deleteIntermediateFiles;
  private final Mode m_mode;


  public static void main(String[] args) {
    Stopwatch stopwatch = Stopwatch.createStarted();

    try {
      CliHelper cliHelper = new CliHelper(MethodHandles.lookup().lookupClass())
          .addVersion("PharmCAT " + CliUtils.getVersion())
          // named allele matcher args
          .addOption("matcher", "matcher", "run named allele matcher")
          .addOption("vcf", "matcher-vcf", "input VCF file for named allele matcher", false, "file")
          .addOption("ma", "matcher-all-results", "return all possible diplotypes, not just top hits")
          .addOption("md", "matcher-definitions-dir", "directory containing named allele definitions (JSON files)", false, "dir")
          .addOption("matcherHtml", "matcher-save-html", "save named allele matcher results as HTML")

          // phenotyper args
          .addOption("phenotyper", "phenotyper", "run phenotyper")
          .addOption("pi", "phenotyper-input", "JSON results from named allele matcher", false, "file")
          .addOption("po", "phenotyper-outside-call-file", "path to an outside call file (TSV)", false, "file")

          // reporter args
          .addOption("reporter", "reporter", "run reporter")
          .addOption("ri", "reporter-input", "JSON results from phenotyper", false, "file")
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
          return;
        }
        if (cliHelper.hasOption("pi")) {
          System.out.println("Cannot specify phenotyper-input (-pi) if running named allele matcher");
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
          return;
        }
      }

      PharmCAT pharmcat = new PharmCAT(config.runMatcher, vcfFile, config.definitionReader,
          config.topCandidateOnly, config.callCyp2d6, config.findCombinations, config.matcherHtml,
          config.runPhenotyper, phenotyperInputFile, phenotyperOutsideCallsFile,
          config.runReporter, reporterInputFile, config.reporterTitle,
          config.reporterSources, config.reporterCompact, config.reporterJson,
          config.outputDir, config.baseFilename, config.deleteIntermediateFiles, Mode.CLI);

      if (!pharmcat.execute()) {
        cliHelper.printHelp();
        return;
      }

      System.out.println("Done.");
      if (cliHelper.isVerbose()) {
        System.out.println("Took " + stopwatch.elapsed(TimeUnit.MILLISECONDS) + "ms");
      }

    } catch (CliHelper.InvalidPathException | ReportableException ex) {
      System.out.println(ex.getMessage());
    } catch (Exception e) {
      e.printStackTrace();
    }
  }


  public PharmCAT(Path vcfFile) throws IOException {
    this(true, vcfFile, null, true, false, false, false,
        true, null, null,
        true, null, null, null, false, false,
        null, null, true, Mode.CLI);
  }

  public PharmCAT(boolean runMatcher, Path vcfFile, @Nullable DefinitionReader definitionReader,
      boolean topCandidateOnly, boolean callCyp2d6, boolean findCombinations, boolean matcherHtml,
      boolean runPhenotyper, @Nullable Path phenotyperInputFile, @Nullable Path phenotyperOutsideCallsFile,
      boolean runReporter, @Nullable Path reporterInputFile, @Nullable String reporterTitle,
      @Nullable List<DataSource> reporterSources, boolean reporterCompact, boolean reporterJson,
      @Nullable Path outputDir, @Nullable String baseFilename, boolean deleteIntermediateFiles, Mode mode)
      throws IOException {
    m_runMatcher = runMatcher;
    if (runMatcher) {
      m_vcfFile = vcfFile;
      m_matcherJsonFile = getOutputFile(m_vcfFile, outputDir, baseFilename, ".match.json");
      if (definitionReader != null) {
        m_definitionReader = definitionReader;
      } else {
        m_definitionReader = new DefinitionReader();
        m_definitionReader.read(DataManager.DEFAULT_DEFINITION_DIR);
      }
      m_topCandidateOnly = topCandidateOnly;
      m_callCyp2d6 = callCyp2d6;
      m_findCombinations = findCombinations;
      if (matcherHtml) {
        m_matcherHtmlFile = getOutputFile(m_vcfFile, outputDir, baseFilename, ".match.html");
      }
    }

    m_runPhenotyper = runPhenotyper;
    if (runPhenotyper) {
      m_phenotyperInputFile = phenotyperInputFile;
      m_phenotyperOutsideCallsFile = phenotyperOutsideCallsFile;
      Path inputFile = m_matcherJsonFile;
      if (m_phenotyperInputFile != null) {
        inputFile = m_phenotyperInputFile;
      } else if (m_phenotyperOutsideCallsFile != null) {
        inputFile = m_phenotyperOutsideCallsFile;
      }
      if (inputFile == null) {
        throw new IllegalStateException("No phenotyper input file");
      }
      m_phenotyperJsonFile = getOutputFile(inputFile, outputDir, baseFilename, ".phenotype.json");
    }

    m_runReporter = runReporter;
    if (runReporter) {
      m_reporterInputFile = reporterInputFile;
      Path inputFile = m_phenotyperJsonFile;
      if (m_reporterInputFile != null) {
        inputFile = m_reporterInputFile;
      }
      if (inputFile == null) {
        throw new IllegalStateException("No reporter input file");
      }
      m_reporterHtmlFile = getOutputFile(inputFile, outputDir, baseFilename, ".report.html");
      if (reporterJson) {
        m_reporterJsonFile = getOutputFile(inputFile, outputDir, baseFilename, ".report.json");
      }
      m_reporterTitle = reporterTitle;
      m_reporterSources = reporterSources;
      m_reporterCompact = reporterCompact;
    }
    m_deleteIntermediateFiles = deleteIntermediateFiles;
    m_mode = mode;
  }


  /**
   * Run PharmCAT pipeline.
   *
   * @return true if some any PharmCAT module was run, false if nothing was run
   */
  public boolean execute() throws IOException {
    boolean didSomething = false;

    Result matcherResult = null;
    if (m_runMatcher) {
      NamedAlleleMatcher namedAlleleMatcher =
          new NamedAlleleMatcher(m_definitionReader, m_findCombinations, m_topCandidateOnly, m_callCyp2d6)
              .printWarnings();
      matcherResult = namedAlleleMatcher.call(m_vcfFile);

      if (m_mode == Mode.CLI) {
        if (!m_deleteIntermediateFiles) {
          System.out.println("Saving named allele matcher JSON results to " + m_matcherJsonFile);
        }
        if (m_matcherHtmlFile != null) {
          System.out.println("Saving named allele matcher HTML results to " + m_matcherHtmlFile);
        }
      }
      namedAlleleMatcher.saveResults(matcherResult, m_matcherJsonFile, m_matcherHtmlFile);

      didSomething = true;
    }

    if (m_runPhenotyper) {
      List<GeneCall> calls;
      Map<String,Collection<String>> warnings = new HashMap<>();
      if (matcherResult != null) {
        calls = matcherResult.getGeneCalls();
        warnings = matcherResult.getVcfWarnings();
      } else if (m_phenotyperInputFile != null) {
        Result deserializedMatcherResult = new ResultSerializer().fromJson(m_phenotyperInputFile);
        calls = deserializedMatcherResult.getGeneCalls();
        warnings = deserializedMatcherResult.getVcfWarnings();
      } else {
        calls = new ArrayList<>();
      }

      List<OutsideCall> outsideCalls = new ArrayList<>();
      if (m_phenotyperOutsideCallsFile != null) {
        outsideCalls = OutsideCallParser.parse(m_phenotyperOutsideCallsFile);
      }

      Phenotyper phenotyper = new Phenotyper(calls, outsideCalls, warnings);
      if (m_mode == Mode.CLI) {
        System.out.println("Saving phenotyper JSON results to " + m_phenotyperJsonFile);
      }
      phenotyper.write(m_phenotyperJsonFile);
      didSomething = true;
    }

    if (m_runReporter) {
      Path inputFile = m_phenotyperJsonFile != null ? m_phenotyperJsonFile : m_reporterInputFile;

      Phenotyper phenotyper = Phenotyper.read(inputFile);
      m_reportContext = new ReportContext(phenotyper.getGeneReports(), m_reporterTitle);
      if (matcherResult != null && matcherResult.getMetadata().isCallCyp2d6()) {
        m_reportContext.addMessage(MessageAnnotation.loadMessage(MessageAnnotation.MSG_CYP2D6_MODE));
      }
      if (m_mode == Mode.CLI) {
        if (!m_deleteIntermediateFiles) {
          System.out.println("Saving reporter HTML results to " + m_reporterHtmlFile);
        }
        if (m_reporterJsonFile != null) {
          System.out.println("Saving reporter JSON results to " + m_reporterJsonFile);
        }
      }
      new HtmlFormat(m_reporterHtmlFile, m_mode == Mode.TEST)
          .sources(m_reporterSources)
          .compact(m_reporterCompact)
          .write(m_reportContext);
      if (m_reporterJsonFile != null) {
        new JsonFormat(m_reporterJsonFile)
            .write(m_reportContext);
      }
      didSomething = true;
    }


    if (m_deleteIntermediateFiles) {
      Files.deleteIfExists(m_matcherJsonFile);
      Files.deleteIfExists(m_phenotyperJsonFile);
    }
    return didSomething;
  }


  /**
   * Gets {@link Path} to an output file based on command line arguments.
   */
  private Path getOutputFile(Path inputFile, @Nullable Path outputDir, @Nullable String baseFilename,
      String defaultSuffix) {
    Path dir;
    if (outputDir != null) {
      dir = outputDir;
    } else {
      dir = inputFile.getParent();
    }
    String filename = Objects.requireNonNullElseGet(baseFilename, () -> getBaseFilename(inputFile));
    return dir.resolve(filename + defaultSuffix);
  }

  public static String getBaseFilename(Path inputFile) {
    String filename = FilenameUtils.getBaseName(inputFile.getFileName().toString());
    if (filename.endsWith(".preprocesed")) {
      filename = filename.substring(0, filename.length() - ".preprocesed".length());
    }
    if (filename.endsWith(".match")) {
      filename = filename.substring(0, filename.length() - ".match".length());
    }
    if (filename.endsWith(".outside")) {
      filename = filename.substring(0, filename.length() - ".outside".length());
    }
    if (filename.endsWith(".phenotype")) {
      filename = filename.substring(0, filename.length() - ".phenotype".length());
    }
    return filename;
  }


  public ReportContext getReportContext() {
    return m_reportContext;
  }




  public static class BaseConfig {
    boolean runMatcher = true;
    DefinitionReader definitionReader;
    boolean topCandidateOnly = true;
    boolean findCombinations;
    boolean callCyp2d6;
    boolean matcherHtml;
    boolean runPhenotyper = true;
    boolean runReporter = true;
    String reporterTitle;
    boolean reporterCompact;
    List<DataSource> reporterSources;
    boolean reporterJson;
    Path outputDir;
    String baseFilename;
    boolean deleteIntermediateFiles;


    BaseConfig(CliHelper cliHelper) throws IOException, ReportableException {
      if (cliHelper.hasOption("matcher") || cliHelper.hasOption("phenotyper") || cliHelper.hasOption("reporter")) {
        runMatcher = cliHelper.hasOption("matcher");
        runPhenotyper = cliHelper.hasOption("phenotyper");
        runReporter = cliHelper.hasOption("reporter");
      }

      if (runMatcher) {
        Path namedAlleleDefinitionDir = DataManager.DEFAULT_DEFINITION_DIR;
        if (cliHelper.hasOption("md")) {
          namedAlleleDefinitionDir = cliHelper.getValidDirectory("md", false);
        }
        definitionReader = new DefinitionReader();
        definitionReader.read(namedAlleleDefinitionDir);
        if (definitionReader.getGenes().size() == 0) {
          throw new ReportableException("Did not find any allele definitions at " + namedAlleleDefinitionDir);
        }

        topCandidateOnly = !cliHelper.hasOption("ma");

        if (cliHelper.hasOption("research")) {
          //noinspection UnstableApiUsage
          List<String> types = sf_commaSplitter.splitToStream(Objects.requireNonNull(cliHelper.getValue("research")))
              .map(String::toLowerCase)
              .collect(Collectors.toList());
          if (types.contains("cyp2d6")) {
            types.remove("cyp2d6");
            System.out.println("WARNING: CYP2D6 RESEARCH MODE ENABLED");
            callCyp2d6 = true;
          }
          if (types.contains("combinations") || types.contains("combination")) {
            types.remove("combinations");
            types.remove("combination");
            System.out.println("WARNING: COMBINATIONS RESEARCH MODE ENABLED");
            findCombinations = true;
          }
          if (types.size() > 0) {
            throw new ReportableException("Unrecognized research option: " + String.join(",", types));
          }
        }
        matcherHtml = cliHelper.hasOption("matcherHtml");
      }

      if (runReporter) {
        reporterTitle = cliHelper.getValue("rt");
        reporterCompact = cliHelper.hasOption("rc");
        reporterJson = cliHelper.hasOption("reporterJson");
        if (cliHelper.hasOption("rs")) {
          reporterSources = new ArrayList<>();
          for (String src : sf_commaSplitter.splitToList(Objects.requireNonNull(cliHelper.getValue("rs")))) {
            try {
              reporterSources.add(DataSource.valueOf(src.toUpperCase()));
            } catch (IllegalArgumentException ex) {
              throw new ReportableException("Unknown source: " + src);
            }
          }
        }
      }

      outputDir = null;
      if (cliHelper.hasOption("o")) {
        outputDir = cliHelper.getValidDirectory("o", true);
      }
      baseFilename = cliHelper.getValue("bf");
      deleteIntermediateFiles = cliHelper.hasOption("del");
    }
  }
}

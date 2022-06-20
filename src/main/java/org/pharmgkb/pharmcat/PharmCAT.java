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
  private DefinitionReader m_definitionReader;
  private boolean m_runMatcher = true;
  private Path m_vcfFile;
  private Path m_namedAlleleDefinitionDir = DataManager.DEFAULT_DEFINITION_DIR;
  private boolean m_topCandidateOnly = true;
  private boolean m_findCombinations;
  private boolean m_callCyp2d6;
  private Path m_matcherJsonFile;
  private Path m_matcherHtmlFile;
  private boolean m_matcherHtml = false;

  private boolean m_runPhenotyper = true;
  private Path m_phenotyperInputFile;
  private Path m_phenotyperOutsideCallsFile;
  private Path m_phenotyperJsonFile;

  private boolean m_runReporter = true;
  private Path m_reporterInputFile;
  private String m_reporterTitle;
  private Path m_reporterJsonFile;
  private Path m_reporterHtmlFile;
  private ReportContext m_reportContext;

  private boolean m_deleteIntermediateFiles;
  private boolean m_cliMode;
  private boolean m_testMode;


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

      if (!new PharmCAT(cliHelper).execute()) {
        cliHelper.printHelp();
      } else {
        System.out.println("Done.");
        if (cliHelper.isVerbose()) {
          System.out.println("Took " + stopwatch.elapsed(TimeUnit.MILLISECONDS) + "ms");
        }
      }

    } catch (CliHelper.InvalidPathException | ReportableException ex) {
      System.out.println(ex.getMessage());
    } catch (Exception e) {
      e.printStackTrace();
    }
  }


  /**
   * Constructor for running entire PharmCAT pipeline.
   *
   * @param isTest this value is used in report formats
   */
  public PharmCAT(boolean isTest) {
    m_testMode = isTest;
  }

  /**
   * Constructor for running PharmCAT via command line.
   */
  public PharmCAT(CliHelper cliHelper) throws Exception {
    if (cliHelper.hasOption("matcher") || cliHelper.hasOption("phenotyper") || cliHelper.hasOption("reporter")) {
      m_runMatcher = cliHelper.hasOption("matcher");
      m_runPhenotyper = cliHelper.hasOption("phenotyper");
      m_runReporter = cliHelper.hasOption("reporter");
    }

    if (m_runMatcher) {
      if (cliHelper.hasOption("vcf")) {
        m_vcfFile = cliHelper.getValidFile("vcf", true);
        m_matcherJsonFile = getOutputFile(cliHelper, m_vcfFile, ".match.json");
      } else {
        throw new ReportableException(
            "No input for Named Allele Matcher!",
            "",
            "Please specify a VCF file (-vcf)"
        );
      }
      if (cliHelper.hasOption("ma")) {
        m_topCandidateOnly = !cliHelper.hasOption("ma");
      }
      if (cliHelper.hasOption("md")) {
        m_namedAlleleDefinitionDir = cliHelper.getValidDirectory("md", false);
      }
      if (cliHelper.hasOption("research")) {
        //noinspection UnstableApiUsage
        List<String> types = sf_commaSplitter.splitToStream(Objects.requireNonNull(cliHelper.getValue("research")))
            .map(String::toLowerCase)
            .collect(Collectors.toList());
        if (types.contains("cyp2d6")) {
          types.remove("cyp2d6");
          System.out.println("WARNING: CYP2D6 RESEARCH MODE ENABLED");
          m_callCyp2d6 = true;
        }
        if (types.contains("combinations") || types.contains("combination")) {
          types.remove("combinations");
          types.remove("combination");
          System.out.println("WARNING: COMBINATIONS RESEARCH MODE ENABLED");
          m_findCombinations = true;
        }
        if (types.size() > 0) {
          throw new ReportableException("Unrecognized research option: " + String.join(",", types));
        }
      }
      if (cliHelper.hasOption("matcherHtml")) {
        m_matcherHtmlFile = getOutputFile(cliHelper, m_vcfFile, ".match.html");
      }

      if (cliHelper.hasOption("pi")) {
        throw new ReportableException("Cannot specify phenotyper-input (-pi) if running named allele matcher");
      }
    }

    if (m_runPhenotyper) {
      Path inputFile = m_matcherJsonFile;
      if (cliHelper.hasOption("pi")) {
        m_phenotyperInputFile = cliHelper.getValidFile("pi", true);
        inputFile = m_phenotyperInputFile;
      }
      if (cliHelper.hasOption("po")) {
        m_phenotyperOutsideCallsFile = cliHelper.getValidFile("po", true);
        if (inputFile == null) {
          inputFile = m_phenotyperOutsideCallsFile;
        }
      }

      if (inputFile == null) {
        throw new ReportableException(
            "No input for Phenotyper!",
            "",
            "Either:",
            "  1. Run named allele matcher with VCF input, or",
            "  2. Specify phenotyper-input (-pi) and/or phenotyper-outside-call-file (-po)"
        );
      }
      m_phenotyperJsonFile = getOutputFile(cliHelper, inputFile, ".phenotype.json");

      if (cliHelper.hasOption("ri")) {
        throw new ReportableException("Cannot specify reporter-input (-ri) if running phenotyper");
      }
    }

    if (m_runReporter) {
      Path inputFile = m_phenotyperJsonFile;
      if (cliHelper.hasOption("ri")) {
        m_reporterInputFile = cliHelper.getValidFile("ri", true);
        inputFile = m_reporterInputFile;
      }
      m_reporterTitle = cliHelper.getValue("rt");

      if (inputFile == null) {
        throw new ReportableException(
            "No input for Reporter!",
            "",
            "Either:",
            "  1. Run phenotyper, or",
            "  2. Specify reporter-input (-ri)"
        );
      }
      m_reporterHtmlFile = getOutputFile(cliHelper, inputFile, ".report.html");
      if (cliHelper.hasOption("reporterJson")) {
        m_reporterJsonFile = getOutputFile(cliHelper, inputFile, ".report.json");
      }
    }

    m_deleteIntermediateFiles = cliHelper.hasOption("del");
    m_cliMode = true;
  }

  public PharmCAT includeMatcherHtml() {
    m_matcherHtml = true;
    return this;
  }


  public PharmCAT matchTopCandidateOnly(boolean topCandidateOnly) {
    m_topCandidateOnly = topCandidateOnly;
    return this;
  }

  public PharmCAT matchCombinations(boolean matchCombinations) {
    m_findCombinations = matchCombinations;
    return this;
  }

  public PharmCAT matchCyp2d6(boolean matchCyp2d6) {
    m_callCyp2d6 = matchCyp2d6;
    return this;
  }


  public void execute(Path vcfFile, Path outsideCallsFile) throws ReportableException, IOException {
    m_vcfFile = vcfFile;
    m_phenotyperOutsideCallsFile = outsideCallsFile;

    String baseFilename = FilenameUtils.getBaseName(m_vcfFile.getFileName().toString());
    m_matcherJsonFile = m_vcfFile.getParent().resolve(baseFilename + ".match.json");
    if (m_matcherHtml) {
      m_matcherHtmlFile = m_vcfFile.getParent().resolve(baseFilename + ".match.html");
    }
    m_phenotyperJsonFile = m_vcfFile.getParent().resolve(baseFilename + ".phenotype.json");
    m_reporterTitle = baseFilename;
    m_reporterJsonFile = m_vcfFile.getParent().resolve(baseFilename + ".report.json");
    m_reporterHtmlFile = m_vcfFile.getParent().resolve(baseFilename + ".report.html");
    execute();
  }


  /**
   * Run PharmCAT pipeline.
   *
   * @return true if some any PharmCAT module was run, false if nothing was run
   */
  private boolean execute() throws ReportableException, IOException {
    boolean didSomething = false;

    Result matcherResult = null;
    if (m_runMatcher) {
      if (m_definitionReader == null) {
        m_definitionReader = initalizeDefinitionReader();
      }

      NamedAlleleMatcher namedAlleleMatcher =
          new NamedAlleleMatcher(m_definitionReader, m_findCombinations, m_topCandidateOnly, m_callCyp2d6)
              .printWarnings();
      matcherResult = namedAlleleMatcher.call(m_vcfFile);

      if (m_cliMode) {
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
      if (m_cliMode) {
        System.out.println("Saving phenotyper JSON results to " + m_phenotyperJsonFile);
      }
      phenotyper.write(m_phenotyperJsonFile);
      didSomething = true;
    }

    if (m_runReporter) {
      Path inputFile = m_phenotyperJsonFile != null ? m_phenotyperJsonFile : m_reporterInputFile;
      m_reportContext = new ReportContext(Phenotyper.readGeneReports(inputFile));
      if (m_cliMode) {
        if (!m_deleteIntermediateFiles) {
          System.out.println("Saving reporter HTML results to " + m_reporterHtmlFile);
        }
        if (m_reporterJsonFile != null) {
          System.out.println("Saving reporter JSON results to " + m_reporterJsonFile);
        }
      }
      new HtmlFormat(m_reporterHtmlFile, m_reporterTitle, m_testMode).write(m_reportContext);
      if (m_reporterJsonFile != null) {
        new JsonFormat(m_reporterJsonFile, m_reporterTitle).write(m_reportContext);
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
  private Path getOutputFile(CliHelper cliHelper, Path inputFile, String defaultSuffix) throws IOException {
    Path outputDir;
    if (cliHelper.hasOption("o")) {
      outputDir = cliHelper.getValidDirectory("o", true);
    } else {
      outputDir = inputFile.getParent();
    }
    String baseFilename;
    if (cliHelper.hasOption("bf")) {
      baseFilename = cliHelper.getValue("bf");
    } else {
      baseFilename = FilenameUtils.getBaseName(inputFile.getFileName().toString());
      if (baseFilename.endsWith(".match")) {
        baseFilename = baseFilename.substring(0, baseFilename.length() - ".match".length());
      }
      if (baseFilename.endsWith(".phenotype")) {
        baseFilename = baseFilename.substring(0, baseFilename.length() - ".phenotype".length());
      }
    }
    return outputDir.resolve(baseFilename + defaultSuffix);
  }


  private DefinitionReader initalizeDefinitionReader() throws ReportableException, IOException {
    DefinitionReader definitionReader = new DefinitionReader();
    definitionReader.read(m_namedAlleleDefinitionDir);
    if (definitionReader.getGenes().size() == 0) {
      throw new ReportableException("Did not find any allele definitions at " + m_namedAlleleDefinitionDir);
    }
    return definitionReader;
  }

  public ReportContext getReportContext() {
    return m_reportContext;
  }
}

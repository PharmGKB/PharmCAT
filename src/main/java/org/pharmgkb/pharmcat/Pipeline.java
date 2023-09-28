package org.pharmgkb.pharmcat;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.concurrent.Callable;
import org.apache.commons.io.FileUtils;
import org.checkerframework.checker.nullness.qual.Nullable;
import org.pharmgkb.common.util.AnsiConsole;
import org.pharmgkb.pharmcat.haplotype.NamedAlleleMatcher;
import org.pharmgkb.pharmcat.haplotype.ResultSerializer;
import org.pharmgkb.pharmcat.haplotype.model.GeneCall;
import org.pharmgkb.pharmcat.phenotype.OutsideCallParser;
import org.pharmgkb.pharmcat.phenotype.Phenotyper;
import org.pharmgkb.pharmcat.phenotype.model.OutsideCall;
import org.pharmgkb.pharmcat.reporter.ReportContext;
import org.pharmgkb.pharmcat.reporter.format.HtmlFormat;
import org.pharmgkb.pharmcat.reporter.format.JsonFormat;
import org.pharmgkb.pharmcat.reporter.model.DataSource;


/**
 * This class runs the PharmCAT modules in a single pipeline.
 *
 * @author Mark Woon
 */
public class Pipeline implements Callable<PipelineResult> {
  public enum Mode {
    /**
     * Default mode.  Prints informative messages to console.
     */
    CLI,
    BATCH,
    /**
     * In test mode, PharmCAT tries not to include version/timestamp in output to simplify diffing.
     */
    TEST
  }
  private final Env m_env;
  private final boolean m_runMatcher;
  private VcfFile m_vcfFile;
  private String m_sampleId;
  private boolean m_topCandidateOnly = true;
  private boolean m_findCombinations;
  private boolean m_callCyp2d6;
  private Path m_matcherJsonFile;
  private Path m_matcherHtmlFile;
  /** True if VCF file only contains a single sample. */
  private final boolean m_singleSample;

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
  private final boolean m_verbose;
  private Path m_baseDir;
  private String m_basename;
  private String m_displayName;
  private final String m_displayCount;


  public Pipeline(Env env,
      boolean runMatcher, @Nullable VcfFile vcfFile, @Nullable String sampleId, boolean singleSample,
      boolean topCandidateOnly, boolean callCyp2d6, boolean findCombinations, boolean matcherHtml,
      boolean runPhenotyper, @Nullable Path phenotyperInputFile, @Nullable Path phenotyperOutsideCallsFile,
      boolean runReporter, @Nullable Path reporterInputFile, @Nullable String reporterTitle,
      @Nullable List<DataSource> reporterSources, boolean reporterCompact, boolean reporterJson, boolean reporterHtml,
      @Nullable Path outputDir, @Nullable String baseFilename, boolean deleteIntermediateFiles,
      Mode mode, @Nullable String displayCount, boolean verbose) throws ReportableException {
    m_env = env;

    m_runMatcher = runMatcher;
    m_baseDir = outputDir;
    if (runMatcher) {
      m_vcfFile = Objects.requireNonNull(vcfFile);
      m_sampleId = sampleId;
      generateBasename(baseFilename, vcfFile.getFile(), sampleId, singleSample);
      if (m_baseDir == null) {
        m_baseDir = m_vcfFile.getFile().getParent();
      }
      m_matcherJsonFile = m_baseDir.resolve(m_basename + BaseConfig.MATCHER_SUFFIX + ".json");
      m_topCandidateOnly = topCandidateOnly;
      m_callCyp2d6 = callCyp2d6;
      m_findCombinations = findCombinations;
      if (matcherHtml) {
        m_matcherHtmlFile = m_baseDir.resolve(m_basename + BaseConfig.MATCHER_SUFFIX + ".html");
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
      generateBasename(baseFilename, inputFile, sampleId, singleSample);
      if (m_baseDir == null) {
        m_baseDir = inputFile.getParent();
      }
      m_phenotyperJsonFile = m_baseDir.resolve(m_basename + BaseConfig.PHENOTYPER_SUFFIX + ".json");
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
      generateBasename(baseFilename, inputFile, sampleId, singleSample);
      if (m_baseDir == null) {
        m_baseDir = inputFile.getParent();
      }
      if (reporterHtml) {
        m_reporterHtmlFile = m_baseDir.resolve(m_basename + BaseConfig.REPORTER_SUFFIX + ".html");
      }
      if (reporterJson) {
        m_reporterJsonFile = m_baseDir.resolve(m_basename + BaseConfig.REPORTER_SUFFIX + ".json");
      }
      m_reporterTitle = reporterTitle;
      if (m_reporterTitle == null) {
        m_reporterTitle = m_basename;
      }
      m_reporterSources = reporterSources;
      m_reporterCompact = reporterCompact;
    }
    m_deleteIntermediateFiles = deleteIntermediateFiles;
    m_mode = mode;
    m_verbose = verbose;
    m_singleSample = singleSample;
    m_displayCount = displayCount;
    Objects.requireNonNull(m_baseDir);
    Objects.requireNonNull(m_basename);
  }


  public @Nullable String getSampleId() {
    return m_sampleId;
  }

  public String getBasename() {
    return m_basename;
  }

  private void generateBasename(String baseFilename, Path inputFile, String sampleId, boolean singleSample) {
    if (m_baseDir == null) {
      m_baseDir = inputFile.getParent();
    }
    if (m_basename != null) {
      return;
    }
    //noinspection ReplaceNullCheck
    if (baseFilename != null) {
      m_basename = baseFilename;
    } else {
      m_basename = BaseConfig.getBaseFilename(inputFile);
    }
    m_displayName = m_basename;
    if (!singleSample &&  sampleId != null && !m_basename.equals(m_sampleId) && !m_basename.startsWith(sampleId + ".") &&
        !m_basename.contains("." + sampleId + ".")) {
      m_basename += "." + sampleId;
    }
    if (sampleId != null) {
      m_displayName = "sample " + sampleId + " in " + m_displayName;
    }
  }


  /**
   * Run PharmCAT pipeline.
   *
   * @throws IOException if there is an error and not running in {@link Mode#BATCH} mode
   */
  @Override
  public PipelineResult call() throws IOException {
    boolean didSomething = false;
    boolean batchDisplayMode = !m_singleSample || m_mode == Mode.BATCH;

    if (batchDisplayMode) {
      StringBuilder builder = new StringBuilder("+ ");
      if (m_displayCount != null) {
        builder.append(m_displayCount)
            .append(" ");
      }
      builder.append("Starting ")
          .append(m_displayName);
      if (m_verbose) {
        builder.append(" (inputs: ")
            .append(getInputDescription())
            .append(")");
      }
      System.out.println(builder);
    }

    try {
      List<String> output = new ArrayList<>();
      org.pharmgkb.pharmcat.haplotype.model.Result matcherResult = null;
      if (m_runMatcher) {
        NamedAlleleMatcher namedAlleleMatcher = new NamedAlleleMatcher(m_env, m_env.getDefinitionReader(),
            m_findCombinations, m_topCandidateOnly, m_callCyp2d6);
        if (!batchDisplayMode) {
          namedAlleleMatcher.printWarnings();
        }
        matcherResult = namedAlleleMatcher.call(m_vcfFile, m_sampleId);

        if (matcherResult.getVcfWarnings() != null &&
            !matcherResult.getVcfWarnings().isEmpty()) {
          Path txtFile = m_matcherJsonFile.getParent()
              .resolve(m_basename + BaseConfig.MATCHER_SUFFIX + "_warnings.txt");
          try (PrintWriter writer = new PrintWriter(Files.newBufferedWriter(txtFile))) {
            Map<String, Collection<String>> warningsMap = matcherResult.getVcfWarnings();
            warningsMap.keySet()
                .forEach(key -> {
                  writer.println(key);
                  warningsMap.get(key)
                      .forEach(msg -> writer.println("\t" + msg));
                });
          }
          output.add(AnsiConsole.styleWarning("Saving VCF warnings to " + txtFile));
        }

        if (!m_deleteIntermediateFiles || !m_runPhenotyper) {
          if (!batchDisplayMode) {
            output.add("Saving named allele matcher JSON results to " + m_matcherJsonFile);
            if (m_matcherHtmlFile != null) {
              output.add("Saving named allele matcher HTML results to " + m_matcherHtmlFile);
            }
          }
          namedAlleleMatcher.saveResults(matcherResult, m_matcherJsonFile, m_matcherHtmlFile);
        }

        didSomething = true;
      }

      Phenotyper phenotyper = null;
      if (m_runPhenotyper) {
        List<GeneCall> calls;
        Map<String, Collection<String>> warnings = new HashMap<>();
        if (matcherResult != null) {
          calls = matcherResult.getGeneCalls();
          warnings.putAll(matcherResult.getVcfWarnings());
        } else if (m_phenotyperInputFile != null) {
          org.pharmgkb.pharmcat.haplotype.model.Result deserializedMatcherResult = new ResultSerializer()
              .fromJson(m_phenotyperInputFile);
          calls = deserializedMatcherResult.getGeneCalls();
          warnings.putAll(deserializedMatcherResult.getVcfWarnings());
        } else {
          calls = new ArrayList<>();
        }

        List<OutsideCall> outsideCalls = new ArrayList<>();
        if (m_phenotyperOutsideCallsFile != null) {
          for (OutsideCall call : OutsideCallParser.parse(m_phenotyperOutsideCallsFile)) {
            if (!m_env.hasGene(call.getGene())) {
              String msg = "Discarded outside call for " + call.getGene() + " because it is not supported by PharmCAT.";
              output.add(AnsiConsole.styleWarning(msg));
              continue;
            }
            if (!m_env.isActivityScoreGene(call.getGene())) {
              if (call.getDiplotype() == null && call.getPhenotype() == null) {
                String msg = call.getGene() + " is not an activity score gene but has outside call with only an " +
                    "activity score.  PharmCAT will not be able to provide any recommendations based on this gene.";
                output.add(AnsiConsole.styleWarning(msg));
              }
            }
            outsideCalls.add(call);
          }
        }

        phenotyper = new Phenotyper(m_env, calls, outsideCalls, warnings);
        if (!m_deleteIntermediateFiles || !m_runReporter) {
          if (!batchDisplayMode) {
            output.add("Saving phenotyper JSON results to " + m_phenotyperJsonFile);
          }
          phenotyper.write(m_phenotyperJsonFile);
        }
        didSomething = true;
      }

      if (m_runReporter) {
        if (phenotyper == null) {
          Path inputFile = m_phenotyperJsonFile != null ? m_phenotyperJsonFile : m_reporterInputFile;
          phenotyper = Phenotyper.read(inputFile);
        }
        m_reportContext = new ReportContext(m_env, phenotyper.getGeneReports(), m_reporterTitle);
        if (m_reporterHtmlFile != null) {
          if (!batchDisplayMode) {
            output.add("Saving reporter HTML results to " + m_reporterHtmlFile);
          }
          new HtmlFormat(m_reporterHtmlFile, m_env, m_mode == Mode.TEST)
              .sources(m_reporterSources)
              .compact(m_reporterCompact)
              .write(m_reportContext);
        }
        if (m_reporterJsonFile != null) {
          if (!batchDisplayMode) {
            output.add("Saving reporter JSON results to " + m_reporterJsonFile);
          }
          new JsonFormat(m_reporterJsonFile, m_env)
              .write(m_reportContext);
        }
        didSomething = true;
      }


      if (m_deleteIntermediateFiles) {
        if (m_matcherJsonFile != null) {
          Files.deleteIfExists(m_matcherJsonFile);
        }
        if (m_phenotyperJsonFile != null) {
          Files.deleteIfExists(m_phenotyperJsonFile);
        }
      }

      StringBuilder builder = new StringBuilder();
      if (batchDisplayMode) {
        builder.append("- ");
        if (m_displayCount != null) {
          builder.append(m_displayCount)
              .append(" ");
        }
        builder.append("Finished processing ")
            .append(m_displayName);
        if (m_verbose) {
          builder.append(" (current memory usage: ")
              .append(FileUtils.byteCountToDisplaySize(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory()))
              .append(")");
        }
        builder.append(System.lineSeparator());
        output.forEach(o -> builder.append("  * ").append(o).append(System.lineSeparator()));
        System.out.print(builder);
      } else {
        output.forEach(System.out::println);
      }
      return new PipelineResult((didSomething ? PipelineResult.Status.SUCCESS : PipelineResult.Status.NOOP), m_basename,
          m_sampleId);

    } catch (Exception ex) {
      if (!m_singleSample || batchDisplayMode) {
        System.err.println("Error with " + m_displayName + ":");
        //noinspection CallToPrintStackTrace
        ex.printStackTrace();
        Path txtFile = m_baseDir.resolve(m_basename + ".ERROR.txt");
        try (PrintWriter writer = new PrintWriter(Files.newBufferedWriter(txtFile))) {
          ex.printStackTrace(writer);
        }
        return new PipelineResult(PipelineResult.Status.FAILURE, m_basename, m_sampleId);
      }
      throw ex;
    }
  }


  public ReportContext getReportContext() {
    return m_reportContext;
  }


  private String getInputDescription() {
    StringBuilder builder = new StringBuilder();
    if (m_vcfFile != null) {
      builder.append(m_vcfFile.getFile().getFileName());
    }
    if (m_phenotyperInputFile != null) {
      if (!builder.isEmpty()) {
        builder.append(", ");
      }
      builder.append(m_phenotyperInputFile.getFileName());
    }
    if (m_phenotyperOutsideCallsFile != null) {
      if (!builder.isEmpty()) {
        builder.append(", ");
      }
      builder.append(m_phenotyperOutsideCallsFile.getFileName());
    }
    if (m_reporterInputFile != null) {
      if (!builder.isEmpty()) {
        builder.append(", ");
      }
      builder.append(m_reporterInputFile.getFileName());
    }
    return builder.toString();
  }

  @Override
  public String toString() {
    return m_displayName;
  }
}

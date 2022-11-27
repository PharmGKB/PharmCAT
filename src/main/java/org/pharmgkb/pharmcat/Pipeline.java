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
import java.util.concurrent.Callable;
import com.google.common.base.Preconditions;
import org.checkerframework.checker.nullness.qual.Nullable;
import org.pharmgkb.pharmcat.haplotype.NamedAlleleMatcher;
import org.pharmgkb.pharmcat.haplotype.ResultSerializer;
import org.pharmgkb.pharmcat.haplotype.model.GeneCall;
import org.pharmgkb.pharmcat.haplotype.model.Result;
import org.pharmgkb.pharmcat.phenotype.OutsideCallParser;
import org.pharmgkb.pharmcat.phenotype.Phenotyper;
import org.pharmgkb.pharmcat.phenotype.model.OutsideCall;
import org.pharmgkb.pharmcat.reporter.ReportContext;
import org.pharmgkb.pharmcat.reporter.format.HtmlFormat;
import org.pharmgkb.pharmcat.reporter.format.JsonFormat;
import org.pharmgkb.pharmcat.reporter.model.DataSource;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


/**
 * This class runs the PharmCAT modules in a single pipeline.
 *
 * @author Mark Woon
 */
public class Pipeline implements Callable<Boolean> {
  private static final Logger sf_logger = LoggerFactory.getLogger(MethodHandles.lookup().lookupClass());
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
  };
  private final Env m_env;
  private final boolean m_runMatcher;
  private Path m_vcfFile;
  private String m_sampleId;
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


  public Pipeline(Env env,
      boolean runMatcher, Path vcfFile, @Nullable String sampleId,
      boolean topCandidateOnly, boolean callCyp2d6, boolean findCombinations, boolean matcherHtml,
      boolean runPhenotyper, @Nullable Path phenotyperInputFile, @Nullable Path phenotyperOutsideCallsFile,
      boolean runReporter, @Nullable Path reporterInputFile, @Nullable String reporterTitle,
      @Nullable List<DataSource> reporterSources, boolean reporterCompact, boolean reporterJson,
      @Nullable Path outputDir, @Nullable String baseFilename, boolean deleteIntermediateFiles, Mode mode)
      throws ReportableException {
    m_env = env;

    m_runMatcher = runMatcher;
    if (runMatcher) {
      m_vcfFile = vcfFile;
      m_sampleId = sampleId;
      if (baseFilename == null) {
        baseFilename = BaseConfig.getBaseFilename(vcfFile);
      }
      if (sampleId != null && !baseFilename.contains("." + sampleId)) {
        baseFilename += "." + sampleId;
      }
      m_matcherJsonFile = getOutputFile(m_vcfFile, outputDir, baseFilename, ".match.json");
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
      if (baseFilename == null) {
        baseFilename = BaseConfig.getBaseFilename(inputFile);
      }
      if (sampleId != null && !baseFilename.contains("." + sampleId)) {
        baseFilename += "." + sampleId;
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
      if (baseFilename == null) {
        baseFilename = BaseConfig.getBaseFilename(inputFile);
      }
      if (sampleId != null && !baseFilename.contains("." + sampleId)) {
        baseFilename += "." + sampleId;
      }
      m_reporterHtmlFile = getOutputFile(inputFile, outputDir, baseFilename, ".report.html");
      if (reporterJson) {
        m_reporterJsonFile = getOutputFile(inputFile, outputDir, baseFilename, ".report.json");
      }
      m_reporterTitle = reporterTitle;
      if (m_reporterTitle == null) {
        m_reporterTitle = baseFilename;
      }
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
  @Override
  public Boolean call() throws IOException {
    boolean didSomething = false;

    Result matcherResult = null;
    if (m_runMatcher) {
      NamedAlleleMatcher namedAlleleMatcher =
          new NamedAlleleMatcher(m_env.getDefinitionReader(), m_findCombinations, m_topCandidateOnly, m_callCyp2d6);
      if (m_mode == Mode.CLI) {
        namedAlleleMatcher.printWarnings();
      }
      matcherResult = namedAlleleMatcher.call(m_vcfFile, m_sampleId);

      if (m_mode == Mode.CLI) {
        if (!m_deleteIntermediateFiles) {
          System.out.println("Saving named allele matcher JSON results to " + m_matcherJsonFile);
        }
        if (m_matcherHtmlFile != null) {
          System.out.println("Saving named allele matcher HTML results to " + m_matcherHtmlFile);
        }
      }
      if (!m_deleteIntermediateFiles || !m_runPhenotyper) {
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
        for (OutsideCall call : OutsideCallParser.parse(m_phenotyperOutsideCallsFile)) {
          if (!m_env.hasGene(call.getGene())) {
            String msg = "Discarded outside call for " + call.getGene() + " because it is not supported by PharmCAT.";
            warnings.put(call.getGene(), List.of(msg));
            sf_logger.warn(msg);
            continue;
          }
          if (!m_env.isActivityScoreGene(call.getGene())) {
            if (call.getDiplotype() == null && call.getPhenotype() == null) {
              String msg = call.getGene() + " is not an activity score gene but has outside call with only an " +
                  "activity score.  PharmCAT will not be able to provide any recommendations based on this gene.";
              warnings.put(call.getGene(), List.of(msg));
              sf_logger.warn(msg);
            }
          }
          outsideCalls.add(call);
        }
      }

      phenotyper = new Phenotyper(m_env, calls, outsideCalls, warnings);
      if (m_mode == Mode.CLI) {
        System.out.println("Saving phenotyper JSON results to " + m_phenotyperJsonFile);
      }
      if (!m_deleteIntermediateFiles || !m_runReporter) {
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
      if (m_mode == Mode.CLI) {
        if (!m_deleteIntermediateFiles) {
          System.out.println("Saving reporter HTML results to " + m_reporterHtmlFile);
        }
        if (m_reporterJsonFile != null) {
          System.out.println("Saving reporter JSON results to " + m_reporterJsonFile);
        }
      }
      new HtmlFormat(m_reporterHtmlFile, m_env, m_mode == Mode.TEST)
          .sources(m_reporterSources)
          .compact(m_reporterCompact)
          .write(m_reportContext);
      if (m_reporterJsonFile != null) {
        new JsonFormat(m_reporterJsonFile, m_env)
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
  private Path getOutputFile(Path inputFile, @Nullable Path outputDir, String baseFilename, String defaultSuffix) {
    Preconditions.checkNotNull(baseFilename);
    Path dir;
    if (outputDir != null) {
      dir = outputDir;
    } else {
      dir = inputFile.toAbsolutePath().getParent();
    }
    return dir.resolve(baseFilename + defaultSuffix);
  }


  public ReportContext getReportContext() {
    return m_reportContext;
  }

  @Override
  public String toString() {
    if (m_sampleId != null) {
      return m_sampleId;
    }
    if (m_vcfFile != null) {
      return BaseConfig.getBaseFilename(m_vcfFile);
    }
    if (m_phenotyperInputFile != null) {
      return BaseConfig.getBaseFilename(m_phenotyperInputFile);
    }
    if (m_phenotyperOutsideCallsFile != null) {
      return BaseConfig.getBaseFilename(m_phenotyperOutsideCallsFile);
    }
    if (m_reporterInputFile != null) {
      return BaseConfig.getBaseFilename(m_reporterInputFile);
    }
    return super.toString();
  }
}

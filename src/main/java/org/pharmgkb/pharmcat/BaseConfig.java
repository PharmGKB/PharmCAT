package org.pharmgkb.pharmcat;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Objects;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import com.google.common.base.Splitter;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang3.StringUtils;
import org.pharmgkb.common.util.CliHelper;
import org.pharmgkb.pharmcat.reporter.model.DataSource;
import org.pharmgkb.pharmcat.reporter.model.PrescribingGuidanceSource;


/**
 * This class handles command-line argument parsing.
 *
 * @author Mark Woon
 */
public class BaseConfig {
  public static final String VCF_PREPROCESSED_SUFFIX = ".preprocessed";
  public static final String VCF_MISSING_PGX_VAR_SUFFIX = ".missing_pgx_var";
  public static final String MATCHER_SUFFIX = ".match";
  public static final String PHENOTYPER_SUFFIX = ".phenotype";
  public static final String REPORTER_SUFFIX = ".report";
  public static final Pattern OUTSIDE_SUFFIX_PATTERN = Pattern.compile("^(.+)\\.outside\\d*$");
  public static final Pattern OUTSIDE_FILENAME_PATTERN = Pattern.compile("^(.+)\\.outside\\d*\\.tsv$");
  public static final Pattern OUTSIDE_EXTENSION_PATTERN = Pattern.compile("^\\.outside\\d*\\.tsv$");
  private static final Splitter sf_commaSplitter = Splitter.on(",").trimResults().omitEmptyStrings();
  boolean runMatcher = true;
  Path definitionDir;
  boolean topCandidateOnly = true;
  boolean findCombinations;
  boolean callCyp2d6;
  boolean matcherHtml;
  boolean runPhenotyper = true;
  boolean runReporter = true;
  String reporterTitle;
  boolean reporterCompact = true;
  List<PrescribingGuidanceSource> reporterSources;
  boolean reporterJson;
  boolean reporterHtml = true;
  boolean reporterCallsOnlyTsv = false;
  Path outputDir;
  String baseFilename;
  boolean deleteIntermediateFiles;
  boolean verbose;
  SortedSet<String> samples = new TreeSet<>();
  Path sampleMetadataFile;


  BaseConfig(CliHelper cliHelper) throws IOException, ReportableException {
    if (cliHelper.hasOption("matcher") || cliHelper.hasOption("phenotyper") || cliHelper.hasOption("reporter")) {
      runMatcher = cliHelper.hasOption("matcher");
      runPhenotyper = cliHelper.hasOption("phenotyper");
      runReporter = cliHelper.hasOption("reporter");
    }
    if (runMatcher && !runPhenotyper && runReporter) {
      throw new ReportableException("Cannot run matcher and reporter without also running phenotyper.");
    }

    if (cliHelper.hasOption("def")) {
      definitionDir = cliHelper.getValidDirectory("def", false);
    }
    if (cliHelper.hasOption("s") && cliHelper.hasOption("S")) {
      throw new ReportableException("Cannot specify both -s and -S");
    }
    if (cliHelper.hasOption("s")) {
      List<String> values = cliHelper.getValues("s");
      for (String v : values) {
        samples.addAll(sf_commaSplitter.splitToList(v));
      }
    }
    if (cliHelper.hasOption("S")) {
      Path sampleFile = cliHelper.getValidFile("S", true);
      try (BufferedReader reader = Files.newBufferedReader(sampleFile)) {
        String line;
        while ((line = reader.readLine()) != null) {
          String sample = StringUtils.stripToNull(line);
          if (sample == null) {
            continue;
          }
          if (sample.contains(",")) {
            throw new ReportableException("Error: Please remove comma ',' from sample names");
          }
          samples.add(sample);
        }
      }
    }
    if (cliHelper.hasOption("sm")) {
      sampleMetadataFile = cliHelper.getValidFile("sm", true);
    }

    boolean researchMode = false;
    if (runMatcher) {
      topCandidateOnly = !cliHelper.hasOption("ma");

      if (cliHelper.hasOption("research")) {
        List<String> types = sf_commaSplitter.splitToStream(Objects.requireNonNull(cliHelper.getValue("research")))
            .map(String::toLowerCase)
            .collect(Collectors.toList());
        if (types.contains("cyp2d6")) {
          types.remove("cyp2d6");
          System.out.println("WARNING: CYP2D6 RESEARCH MODE ENABLED");
          callCyp2d6 = true;
          researchMode = true;
        }
        if (types.contains("combinations") || types.contains("combination")) {
          types.remove("combinations");
          types.remove("combination");
          System.out.println("WARNING: COMBINATIONS RESEARCH MODE ENABLED");
          findCombinations = true;
          researchMode = true;
        }
        if (!types.isEmpty()) {
          throw new ReportableException("Unrecognized research option: " + String.join(",", types));
        }
      }
      matcherHtml = cliHelper.hasOption("matcherHtml");
    }

    if (runReporter) {
      reporterTitle = cliHelper.getValue("rt");
      reporterCompact = !cliHelper.hasOption("re");
      reporterJson = cliHelper.hasOption("reporterJson");
      reporterCallsOnlyTsv = cliHelper.hasOption("reporterCallsOnlyTsv");
      // by default, generate the HTML report (which preserves backwards compatibility)
      // only check for reporterHtml flag if specifying other reporter format outputs
      if (reporterJson || reporterCallsOnlyTsv) {
        reporterHtml = cliHelper.hasOption("reporterHtml");
      }

      if (researchMode) {
        System.out.println("WARNING: FULL REPORTER OUTPUT NOT AVAILABLE IN RESEARCH MODE");
        if (!reporterCallsOnlyTsv) {
          reporterCallsOnlyTsv = true;
        }
        reporterHtml = false;
        reporterJson = false;
      }

      if (cliHelper.hasOption("rs")) {
        reporterSources = new ArrayList<>();
        for (String src : sf_commaSplitter.splitToList(Objects.requireNonNull(cliHelper.getValue("rs")))) {
          try {
            DataSource ds = DataSource.valueOf(src);
            switch (ds) {
              case CPIC:
                reporterSources.add(PrescribingGuidanceSource.CPIC_GUIDELINE);
                break;
              case DPWG:
                reporterSources.add(PrescribingGuidanceSource.DPWG_GUIDELINE);
                break;
              case FDA:
                reporterSources.add(PrescribingGuidanceSource.FDA_LABEL);
                reporterSources.add(PrescribingGuidanceSource.FDA_ASSOC);
                break;
              default:
                throw new ReportableException("Unsupported source: " + src);
            }
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
    verbose = cliHelper.isVerbose();
  }


  public boolean runSample(String sample) {
    if (samples.isEmpty()) {
      return true;
    }
    return samples.contains(sample);
  }

  public static String getBaseFilename(Path inputFile) {
    String filename = FilenameUtils.getBaseName(inputFile.getFileName().toString());
    if (filename.endsWith(".vcf")) {
      // because .vcf might come in as .vcf.bgz or .vcf.gz
      filename = FilenameUtils.getBaseName(filename);
    }
    if (filename.endsWith(VCF_PREPROCESSED_SUFFIX)) {
      filename = filename.substring(0, filename.length() - VCF_PREPROCESSED_SUFFIX.length());
    }
    if (filename.endsWith(MATCHER_SUFFIX)) {
      filename = filename.substring(0, filename.length() - MATCHER_SUFFIX.length());
    }
    Matcher m = BaseConfig.OUTSIDE_SUFFIX_PATTERN.matcher(filename);
    if (m.matches()) {
      filename = m.group(1);
    }
    if (filename.endsWith(PHENOTYPER_SUFFIX)) {
      filename = filename.substring(0, filename.length() - PHENOTYPER_SUFFIX.length());
    }
    if (filename.endsWith(REPORTER_SUFFIX)) {
      filename = filename.substring(0, filename.length() - REPORTER_SUFFIX.length());
    }
    return filename;
  }
}

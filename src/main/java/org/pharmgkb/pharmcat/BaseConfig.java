package org.pharmgkb.pharmcat;

import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Objects;
import java.util.stream.Collectors;
import com.google.common.base.Splitter;
import org.apache.commons.io.FilenameUtils;
import org.pharmgkb.common.util.CliHelper;
import org.pharmgkb.pharmcat.reporter.model.DataSource;


/**
 * This class handles command-line argument parsing.
 *
 * @author Mark Woon
 */
class BaseConfig {
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
      if (cliHelper.hasOption("md")) {
        definitionDir = cliHelper.getValidDirectory("md", false);
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
      reporterCompact = !cliHelper.hasOption("re");
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


  public static String getBaseFilename(Path inputFile) {
    String filename = FilenameUtils.getBaseName(inputFile.getFileName().toString());
    if (filename.endsWith(".preprocessed")) {
      filename = filename.substring(0, filename.length() - ".preprocessed".length());
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
}

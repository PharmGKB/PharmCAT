package org.pharmgkb.pharmcat;

import java.io.IOException;
import java.lang.invoke.MethodHandles;
import java.nio.file.Path;
import java.util.Date;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import javax.annotation.Nonnull;
import javax.annotation.Nullable;
import com.google.common.base.Preconditions;
import org.pharmgkb.common.io.util.CliHelper;
import org.pharmgkb.common.util.PathUtils;
import org.pharmgkb.pharmcat.haplotype.DefinitionReader;
import org.pharmgkb.pharmcat.haplotype.NamedAlleleMatcher;
import org.pharmgkb.pharmcat.haplotype.ResultSerializer;
import org.pharmgkb.pharmcat.haplotype.model.Result;
import org.pharmgkb.pharmcat.reporter.Reporter;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


/**
 * Class to run the PharmCAT tool from input VCF file to final output report.
 *
 * @author Ryan Whaley
 */
public class PharmCAT {
  private static final Logger sf_logger = LoggerFactory.getLogger(MethodHandles.lookup().lookupClass());
  private static final Pattern sf_inputNamePattern = Pattern.compile("(.*)\\.vcf");

  private NamedAlleleMatcher m_namedAlleleMatcher;
  private Reporter m_reporter;
  private Path m_outputDir;

  public static void main(String[] args) {
    CliHelper cliHelper = new CliHelper(MethodHandles.lookup().lookupClass())
        .addOption("c", "call-file", "input call file (VCF)", true, "c")
        .addOption("o", "output-dir", "directory to output to", true, "o")
        .addOption("f", "output-file", "the output name used for ouput file names (will add file extensions), will default to same value as call-file if not specified", false, "f")
        .addOption("n", "annotation-dir", "directory of guideline annotations (JSON files)", true, "n")
        .addOption("l", "alleles-dir", "directory of named allele definitions (JSON files)", false, "l")
        .addOption("a", "astrolabe-file", "path to astrolabe result file (TSV)", false, "a");

    try {
      if (!cliHelper.parse(args)) {
        System.exit(1);
      }

      Path inputFile = cliHelper.getValidFile("c", true);
      Path outputDir = cliHelper.getValidDirectory("o", true);
      Path annoDir = cliHelper.getValidDirectory("n", false);

      Path astrolabeFile = null;
      if (cliHelper.hasOption("a")) {
        astrolabeFile = cliHelper.getPath("a");
      }

      Path allelesDir;
      if (cliHelper.hasOption("l")) {
        allelesDir = cliHelper.getValidDirectory("l", false);
      }
      else {
        allelesDir = PathUtils.getPathToResource("org/pharmgkb/pharmcat/definition/alleles");
      }

      String outputFile = null;
      if (cliHelper.hasOption("f")) {
        outputFile = cliHelper.getValue("f");
      }

      PharmCAT pharmcat = new PharmCAT(outputDir, allelesDir, annoDir);
      pharmcat.execute(inputFile, astrolabeFile, outputFile);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  /**
   * public constructor.
   *
   * Sets up all the necessary supporting objects in order to run the matcher and the reporter.
   *
   * @param outputDir Path to the directory to write output to
   * @param allelesDir Path to the directory where allele definitions are
   * @param annoDir Path to the directory where guideline annotations are
   * @throws IOException can be throwsn if filesystem objects not in proper state
   */
  private PharmCAT(@Nonnull Path outputDir, @Nonnull Path allelesDir, @Nonnull Path annoDir) throws IOException {

    boolean madeDir = outputDir.toFile().mkdirs();
    if (madeDir) {
      sf_logger.info("Directory didn't exist so created "+outputDir);
    }

    Preconditions.checkArgument(outputDir.toFile().exists(), "output directory doesn't exist");
    Preconditions.checkArgument(outputDir.toFile().isDirectory(), "output path is not a directory");

    Preconditions.checkArgument(allelesDir.toFile().exists(), "path to allele definitions does not exist");
    Preconditions.checkArgument(allelesDir.toFile().isDirectory(), "alleles directory path is not a directory");

    Preconditions.checkArgument(annoDir.toFile().exists(), "path to annotations definitions does not exist");
    Preconditions.checkArgument(annoDir.toFile().isDirectory(), "annotations directory path is not a directory");

    DefinitionReader definitionReader = new DefinitionReader();
    definitionReader.read(allelesDir);

    m_namedAlleleMatcher = new NamedAlleleMatcher(definitionReader);
    m_reporter = new Reporter(annoDir);
    m_outputDir = outputDir;

    sf_logger.info("Using alleles: " + allelesDir);
    sf_logger.info("Using annotations: " + annoDir);
    sf_logger.info("Writing to: " + outputDir);
  }

  /**
   * Executes the {@link NamedAlleleMatcher} then the {@link Reporter} on the given sample data
   * @param inputFile the input sample VCF file
   * @param astrolabeFile the optional input astrolabe TSV file
   * @param outputFile the optional name to write the output to
   * @throws Exception can occur from file I/O or unexpected state
   */
  private void execute(@Nonnull Path inputFile, @Nullable Path astrolabeFile, @Nullable String outputFile) throws Exception {
    Preconditions.checkArgument(inputFile.toFile().exists(), "input file does not exist");
    Preconditions.checkArgument(inputFile.toFile().isFile(), "input path is not a file");

    sf_logger.info("Run time: " + new Date());
    String fileRoot = makeFileRoot(inputFile, outputFile);

    Path callPath = m_outputDir.resolve(fileRoot + ".call.json");
    Path reportPath = m_outputDir.resolve(fileRoot + ".report.html");

    Result result = m_namedAlleleMatcher.call(inputFile);
    ResultSerializer resultSerializer = new ResultSerializer();
    resultSerializer.toJson(result, callPath);

    m_reporter.analyze(callPath, astrolabeFile);
    m_reporter.printHtml(reportPath, fileRoot);

    callPath.toFile().deleteOnExit();

    sf_logger.info("Completed");
  }

  /**
   * Determines what to call the output file depending on user input parameters
   * @param inputFile the input VCF file path
   * @param outputFile the optional output file name to use
   * @return a file name to use without extension
   */
  @Nonnull
  private String makeFileRoot(@Nonnull Path inputFile, @Nullable String outputFile) {
    Matcher m = sf_inputNamePattern.matcher(inputFile.getFileName().toString());
    String fileRoot;
    if (outputFile != null) {
      fileRoot = outputFile;
    }
    else if (m.matches()) {
      fileRoot = m.group(1);
    }
    else {
      fileRoot = inputFile.getFileName().toString();
    }
    return fileRoot;
  }
}

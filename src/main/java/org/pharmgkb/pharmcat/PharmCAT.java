package org.pharmgkb.pharmcat;

import java.io.IOException;
import java.lang.invoke.MethodHandles;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Date;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import javax.annotation.Nonnull;
import javax.annotation.Nullable;
import com.google.common.base.Preconditions;
import org.apache.commons.io.FileUtils;
import org.pharmgkb.common.io.util.CliHelper;
import org.pharmgkb.pharmcat.haplotype.DefinitionReader;
import org.pharmgkb.pharmcat.haplotype.NamedAlleleMatcher;
import org.pharmgkb.pharmcat.haplotype.ResultSerializer;
import org.pharmgkb.pharmcat.haplotype.model.Result;
import org.pharmgkb.pharmcat.reporter.Reporter;
import org.pharmgkb.pharmcat.util.DataManager;
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
  private boolean m_keepMatcherOutput = false;
  private boolean m_writeJsonReport = false;

  public static void main(String[] args) {
    CliHelper cliHelper = new CliHelper(MethodHandles.lookup().lookupClass())
        .addOption("vcf", "sample-file", "input call file (VCF)", true, "vcf")
        .addOption("o", "output-dir", "directory to output to", true, "o")
        .addOption("f", "output-file", "the base name used for ouput file names (will add file extensions), will default to same value as call-file if not specified", false, "f")
        .addOption("a", "astrolabe-file", "path to astrolabe result file (TSV)", false, "a")
        // optional data
        .addOption("g", "guidelines-dir", "directory of guideline annotations (JSON files)", false, "n")
        .addOption("na", "alleles-dir", "directory of named allele definitions (JSON files)", false, "l")
        // controls
        .addOption("k", "keep-matcher-files", "flag to keep the intermediary matcher output files")
        .addOption("j", "write-report-json", "flag to write a JSON file of the data used to populate the final report");

    try {
      if (!cliHelper.parse(args)) {
        System.exit(1);
      }

      Path vcfFile = cliHelper.getValidFile("vcf", true);
      Path outputDir = cliHelper.getValidDirectory("o", true);
      Path astrolabeFile = null;
      if (cliHelper.hasOption("a")) {
        astrolabeFile = cliHelper.getPath("a");
      }

      Path guidelinesDir = null;
      if (cliHelper.hasOption("g")) {
        guidelinesDir = cliHelper.getValidDirectory("g", false);
      }
      Path definitionsDir = null;
      if (cliHelper.hasOption("na")) {
        definitionsDir = cliHelper.getValidDirectory("l", false);
      }

      String outputFile = null;
      if (cliHelper.hasOption("f")) {
        outputFile = cliHelper.getValue("f");
      }

      PharmCAT pharmcat = new PharmCAT(outputDir, definitionsDir, guidelinesDir);
      if (cliHelper.hasOption("k")) {
        pharmcat.keepMatcherOutput();
      }

      pharmcat
          .writeJson(cliHelper.hasOption("j"))
          .execute(vcfFile, astrolabeFile, outputFile);

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
   * @param definitionsDir Path to the directory where allele definitions are, null will use default definitions
   * @param guidelinesDir Path to the directory where guideline annotations are, null will use default annotations
   * @throws IOException can be throwsn if filesystem objects not in proper state
   */
  public PharmCAT(@Nonnull Path outputDir, @Nullable Path definitionsDir, @Nullable Path guidelinesDir)
      throws IOException {

    boolean madeDir = outputDir.toFile().mkdirs();
    if (madeDir) {
      sf_logger.info("Directory didn't exist so created {}", outputDir);
    }
    Preconditions.checkArgument(Files.isDirectory(outputDir), "Not a directory: %s", outputDir);

    if (definitionsDir == null) {
      definitionsDir = DataManager.DEFAULT_DEFINITION_DIR;
    }
    Preconditions.checkArgument(Files.isDirectory(definitionsDir), "Not a directory: %s", definitionsDir);

    if (guidelinesDir == null) {
      guidelinesDir = DataManager.DEFAULT_GUIDELINE_DIR;
    }
    Preconditions.checkArgument(Files.isDirectory(guidelinesDir), "Not a directory: %s", guidelinesDir);

    DefinitionReader definitionReader = new DefinitionReader();
    definitionReader.read(definitionsDir);

    m_namedAlleleMatcher = new NamedAlleleMatcher(definitionReader, true, true)
        .printWarnings();
    m_reporter = new Reporter(guidelinesDir);
    m_outputDir = outputDir;

    sf_logger.info("Using alleles: {}", definitionsDir);
    sf_logger.info("Using annotations: {}", guidelinesDir);
    sf_logger.info("Writing to: {}", outputDir);

    System.out.println(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<");
    System.out.println(">>> WARNING: This is pre-release software, not for production use <<<");
    System.out.println(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<");
  }

  /**
   * Executes the {@link NamedAlleleMatcher} then the {@link Reporter} on the given sample data
   * @param vcfFile the input sample VCF file
   * @param astrolabeFile the optional input astrolabe TSV file
   * @param outputFile the optional name to write the output to
   * @throws Exception can occur from file I/O or unexpected state
   */
  public void execute(@Nonnull Path vcfFile, @Nullable Path astrolabeFile, @Nullable String outputFile) throws Exception {
    Preconditions.checkArgument(Files.isRegularFile(vcfFile), "Not a file: %s", vcfFile);

    sf_logger.info("Run time: " + new Date());
    String fileRoot = makeFileRoot(vcfFile, outputFile);

    Path callFile = m_outputDir.resolve(fileRoot + ".call.json");
    if (!m_keepMatcherOutput) {
      callFile.toFile().deleteOnExit();
    }

    Result result = m_namedAlleleMatcher.call(vcfFile);
    ResultSerializer resultSerializer = new ResultSerializer();
    resultSerializer.toJson(result, callFile);
    if (m_keepMatcherOutput) {
      resultSerializer.toHtml(result, m_outputDir.resolve(fileRoot + ".matcher.html"));
    }

    m_reporter.analyze(callFile, astrolabeFile);

    Path reportPath = m_outputDir.resolve(fileRoot + ".report.html");
    Path jsonPath = m_writeJsonReport ? m_outputDir.resolve(fileRoot + ".report.json") : null;
    m_reporter.printHtml(reportPath, fileRoot, jsonPath);

    if (!m_keepMatcherOutput) {
      FileUtils.deleteQuietly(callFile.toFile());
    }

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

  public PharmCAT keepMatcherOutput() {
    m_keepMatcherOutput = true;
    return this;
  }

  /**
   * Determine whether to write reporter JSON output or not
   * @param doWrite true to create a <code>.report.json</code> file as output
   */
  private PharmCAT writeJson(boolean doWrite) {
    m_writeJsonReport = doWrite;
    return this;
  }

  /**
   * Getter for the Reporter class to use for testing
   * @return the current Reporter instance
   */
  public Reporter getReporter() {
    return m_reporter;
  }
}

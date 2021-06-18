package org.pharmgkb.pharmcat;

import java.io.IOException;
import java.lang.invoke.MethodHandles;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import javax.annotation.Nonnull;
import com.google.common.base.Preconditions;
import org.apache.commons.io.FileUtils;
import org.checkerframework.checker.nullness.qual.Nullable;
import org.pharmgkb.common.io.util.CliHelper;
import org.pharmgkb.pharmcat.haplotype.DefinitionReader;
import org.pharmgkb.pharmcat.haplotype.NamedAlleleMatcher;
import org.pharmgkb.pharmcat.haplotype.ResultSerializer;
import org.pharmgkb.pharmcat.haplotype.model.Result;
import org.pharmgkb.pharmcat.phenotype.Phenotyper;
import org.pharmgkb.pharmcat.reporter.Reporter;
import org.pharmgkb.pharmcat.reporter.io.OutsideCallParser;
import org.pharmgkb.pharmcat.reporter.model.OutsideCall;
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

  private final Reporter f_reporter;
  private final Path f_definitionsDir;
  private Path m_outputDir;
  private boolean m_keepMatcherOutput = false;
  private boolean m_writeJsonReport = false;
  private boolean m_writeJsonPheno = false;
  private boolean m_showAllMatches = false;

  public static void main(String[] args) {
    CliHelper cliHelper = new CliHelper(MethodHandles.lookup().lookupClass())
        .addOption("vcf", "sample-file", "input call file (VCF)", true, "vcf")
        .addOption("o", "output-dir", "directory to output to", true, "o")
        .addOption("f", "output-file", "the base name used for ouput file names (will add file extensions), will default to same value as call-file if not specified", false, "f")
        .addOption("a", "outside-call-file", "path to an outside call file (TSV)", false, "a")
        // optional data
        .addOption("na", "alleles-dir", "directory of named allele definitions (JSON files)", false, "l")
        // controls
        .addOption("m", "show-all-matches", "show all matching diplotype pairs, not just top-scoring match")
        .addOption("k", "keep-matcher-files", "flag to keep the intermediary matcher output files")
        .addOption("j", "write-reporter-json", "flag to write a JSON file of the data used to populate the final report")
        .addOption("pj", "write-phenotyper-json", "flag to write a JSON file of the data used in the phenotyper");

    try {
      if (!cliHelper.parse(args)) {
        System.exit(1);
      }

      Path vcfFile = cliHelper.getValidFile("vcf", true);
      Path outputDir = cliHelper.getValidDirectory("o", true);
      Path outsideCallPath = null;
      if (cliHelper.hasOption("a")) {
        outsideCallPath = cliHelper.getPath("a");
      }

      Path definitionsDir = null;
      if (cliHelper.hasOption("na")) {
        definitionsDir = cliHelper.getValidDirectory("l", false);
      }

      String outputFile = null;
      if (cliHelper.hasOption("f")) {
        outputFile = cliHelper.getValue("f");
      }

      PharmCAT pharmcat = new PharmCAT(outputDir, definitionsDir);
      if (cliHelper.hasOption("k")) {
        pharmcat.keepMatcherOutput();
      }
      if (cliHelper.hasOption("m")) {
        pharmcat.showAllMatches();
      }

      pharmcat
          .writeJson(cliHelper.hasOption("j"))
          .writePhenotyperJson(cliHelper.hasOption("pj"))
          .execute(vcfFile, outsideCallPath, outputFile);

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
   * @throws IOException can be throwsn if filesystem objects not in proper state
   */
  public PharmCAT(Path outputDir, @Nullable Path definitionsDir)
      throws IOException {

    boolean madeDir = outputDir.toFile().mkdirs();
    if (madeDir) {
      sf_logger.info("Directory didn't exist so created {}", outputDir);
    }
    Preconditions.checkArgument(Files.isDirectory(outputDir), "Not a directory: %s", outputDir);

    if (definitionsDir == null) {
      f_definitionsDir = DataManager.DEFAULT_DEFINITION_DIR;
    } else {
      f_definitionsDir = definitionsDir;
    }
    Preconditions.checkArgument(Files.isDirectory(f_definitionsDir), "Not a directory: %s", definitionsDir);

    f_reporter = new Reporter();
    setOutputDir(outputDir);

    sf_logger.info("Using alleles: {}", definitionsDir);
    sf_logger.info("Writing to: {}", outputDir);

    System.out.println(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<");
    System.out.println(">>> WARNING: This is pre-release software, not for production use <<<");
    System.out.println(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<");
  }

  /**
   * Executes the {@link NamedAlleleMatcher} then the {@link Reporter} on the given sample data
   * @param vcfFile the input sample VCF file
   * @param outsideCallPath the optional input outside call TSV file
   * @param outputFile the optional name to write the output to
   * @throws Exception can occur from file I/O or unexpected state
   */
  public void execute(Path vcfFile, @Nullable Path outsideCallPath, @Nullable String outputFile) throws Exception {
    Preconditions.checkArgument(Files.isRegularFile(vcfFile), "Not a file: %s", vcfFile);

    sf_logger.info("Run time: " + new Date());
    String fileRoot = makeFileRoot(vcfFile, outputFile);

    Path callFile = m_outputDir.resolve(fileRoot + ".matcher.json");
    if (!m_keepMatcherOutput) {
      callFile.toFile().deleteOnExit();
    }

    DefinitionReader definitionReader = new DefinitionReader();
    definitionReader.read(f_definitionsDir);
    NamedAlleleMatcher namedAlleleMatcher = new NamedAlleleMatcher(definitionReader, true, !m_showAllMatches)
        .printWarnings();
    Result result = namedAlleleMatcher.call(vcfFile);
    ResultSerializer resultSerializer = new ResultSerializer();
    resultSerializer.toJson(result, callFile);
    if (m_keepMatcherOutput) {
      resultSerializer.toHtml(result, m_outputDir.resolve(fileRoot + ".matcher.html"));
    }

    //Load the outside calls if it's available
    List<OutsideCall> outsideCalls = new ArrayList<>();
    if (outsideCallPath != null) {
      Preconditions.checkArgument(Files.exists(outsideCallPath));
      Preconditions.checkArgument(Files.isRegularFile(outsideCallPath));
      outsideCalls = OutsideCallParser.parse(outsideCallPath);
    }

    Phenotyper phenotyper = new Phenotyper(result.getGeneCalls(), outsideCalls);

    f_reporter.analyze(phenotyper.getGeneReports());

    Path reportPath = m_outputDir.resolve(fileRoot + ".report.html");
    Path jsonPath = m_writeJsonReport ? m_outputDir.resolve(fileRoot + ".report.json") : null;
    f_reporter.printHtml(reportPath, fileRoot, jsonPath);

    if (m_writeJsonPheno) {
      phenotyper.write(m_outputDir.resolve(fileRoot + ".phenotyper.json"));
    }


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
  private String makeFileRoot(Path inputFile, @Nullable String outputFile) {
    String fileRoot;
    if (outputFile != null) {
      fileRoot = outputFile;
    } else {
      Matcher m = sf_inputNamePattern.matcher(inputFile.getFileName().toString());
      if (m.matches()) {
        fileRoot = m.group(1);
      }
      else {
        fileRoot = inputFile.getFileName().toString();
      }
    }
    return fileRoot;
  }

  public PharmCAT keepMatcherOutput() {
    m_keepMatcherOutput = true;
    return this;
  }

  public PharmCAT showAllMatches() {
    m_showAllMatches = true;
    return this;
  }

  /**
   * Determine whether to write reporter JSON output or not
   * @param doWrite true to create a <code>.report.json</code> file as output
   */
  public PharmCAT writeJson(boolean doWrite) {
    m_writeJsonReport = doWrite;
    return this;
  }

  public PharmCAT writePhenotyperJson(boolean doWrite) {
    m_writeJsonPheno = doWrite;
    return this;
  }

  /**
   * Getter for the Reporter class to use for testing
   * @return the current Reporter instance
   */
  public Reporter getReporter() {
    return f_reporter;
  }

  /**
   * Sets the output directory and does some sanity checking
   * @param outputDir a nonnull path to an existing directory
   */
  public void setOutputDir(Path outputDir) {
    Preconditions.checkNotNull(outputDir);
    Preconditions.checkArgument(
        outputDir.toFile().exists() || outputDir.toFile().isDirectory(),
        "Specified output isn't a directory: " + outputDir.toAbsolutePath());
    
    m_outputDir = outputDir;
  }
}

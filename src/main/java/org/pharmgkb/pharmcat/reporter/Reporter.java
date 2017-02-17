package org.pharmgkb.pharmcat.reporter;

import java.io.IOException;
import java.lang.invoke.MethodHandles;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.stream.Collectors;
import javax.annotation.Nonnull;
import javax.annotation.Nullable;
import com.google.common.base.Preconditions;
import org.pharmgkb.common.io.util.CliHelper;
import org.pharmgkb.pharmcat.annotation.AnnotationReader;
import org.pharmgkb.pharmcat.haplotype.model.GeneCall;
import org.pharmgkb.pharmcat.reporter.io.JsonFileLoader;
import org.pharmgkb.pharmcat.reporter.io.MarkdownWriter;
import org.pharmgkb.pharmcat.reporter.model.GuidelinePackage;
import org.pharmgkb.pharmcat.reporter.model.result.GuidelineReport;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


/**
 * This is contains the main class for running the reporting tool. It's responsible for taking input of all the
 * necessary data files, parsing them, and running the reporter codes.
 *
 * This can be run both on the command line and procedurally.
 *
 * @author greytwist
 * @author Ryan Whaley
 */
public class Reporter {
  private static final Logger sf_logger = LoggerFactory.getLogger(MethodHandles.lookup().lookupClass());
  private AnnotationReader m_annotationReader;
  private List<Path> m_annotationFiles = null;
  private ReportContext m_reportContext = null;

  /**
   * main
   * @param args command line args
   */
  public static void main(String[] args) {

    CliHelper cliHelper = new CliHelper(MethodHandles.lookup().lookupClass())
        .addOption("d", "annotations-dir", "directory of allele definition files", true, "d")
        .addOption("c", "call-file", "named allele call file", true, "c")
        .addOption("o", "output-file", "file to write report output to", true, "o")
        ;

    try {
      if (!cliHelper.parse(args)) {
        System.exit(1);
      }

      Path annotationsDir = cliHelper.getValidDirectory("d", false);
      Path callFile = cliHelper.getValidFile("c", true);
      Path outputFile = cliHelper.getPath("o");

      new Reporter(annotationsDir)
          .analyze(callFile)
          .printMarkdown(outputFile);

    } catch (Exception ex) {
      ex.printStackTrace();
    }
  }

  /**
   * public constructor. start a new reporter based on annotation data found in the given <code>annotationsDir</code>.
   *
   * @param annotationsDir directory of annotation files
   */
  public Reporter(@Nonnull Path annotationsDir) throws IOException {

    Preconditions.checkNotNull(annotationsDir);
    Preconditions.checkArgument(Files.exists(annotationsDir));
    Preconditions.checkArgument(Files.isDirectory(annotationsDir));

    m_annotationReader = new AnnotationReader();
    m_annotationReader.read(annotationsDir);

    m_annotationFiles = Files.list(annotationsDir)
        .filter(f -> f.getFileName().toString().endsWith(".json"))
        .collect(Collectors.toList());
    if (m_annotationFiles.size() == 0) {
      throw new IOException("No annotation definitions to read from");
    }
  }

  /**
   * Run the actual report process. Parse the input file, do the matching, and write the report files.
   *
   * @param callFile file of haplotype calls
   */
  public Reporter analyze(@Nonnull Path callFile) throws Exception {
    Preconditions.checkNotNull(callFile);
    Preconditions.checkArgument(Files.exists(callFile));
    Preconditions.checkArgument(Files.isRegularFile(callFile));

    //Generate class used for loading JSON into
    JsonFileLoader loader = new JsonFileLoader();

    //Load the haplotype json, this is pointed at a test json and will likely break when meeting real
    // requiring some if not all rewriting
    List<GeneCall> calls = loader.loadHaplotypeGeneCalls(callFile);

    //Load the gene drug interaction list. This currently only handles single gene drug m_guidelineFiles and will require updating to handle multi gene drug interaction
    List<GuidelinePackage> guidelines = loader.loadGuidelines(m_annotationFiles);

    //This is the primary work flow for generating the report where calls are matched to exceptions and drug gene m_guidelineFiles based on reported haplotypes
    m_reportContext = new ReportContext(calls, guidelines);

    return this;
  }

  /**
   * Print a Markdown file of compiled report data
   * @param reportFile directory to write output to
   */
  public void printMarkdown(@Nonnull Path reportFile) throws IOException {
    Preconditions.checkNotNull(reportFile);

    new MarkdownWriter(reportFile)
        .print(m_reportContext);
  }

  /**
   * Expose the guideline reports for testing purposes
   */
  @Nullable
  protected List<GuidelineReport> getGuidelineReports() {
    if (m_reportContext == null) {
      return null;
    }
    return m_reportContext.getGuidelineResults();
  }
}

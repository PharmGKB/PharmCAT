package org.pharmgkb.pharmcat.reporter;

import java.io.File;
import java.io.IOException;
import java.lang.invoke.MethodHandles;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;
import javax.annotation.Nonnull;
import javax.annotation.Nullable;
import com.google.common.base.Preconditions;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.pharmgkb.pharmcat.haplotype.model.GeneCall;
import org.pharmgkb.pharmcat.reporter.io.JsonFileLoader;
import org.pharmgkb.pharmcat.reporter.io.MarkdownWriter;
import org.pharmgkb.pharmcat.reporter.model.DosingGuideline;
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

  private List<File> m_annotationFiles = null;
  private DataUnifier m_dataUnifier = null;

  /**
   * main
   * @param args command line args
   */
  public static void main(String[] args) throws Exception {
    Options options = new Options();

    options.addOption(new Option("annotationsDir", true, "required - directory holding all the annotations files"));
    options.addOption(new Option("callFile", true, "required - file from the Haplotyper"));
    options.addOption(new Option("reportFile", true, "required - file to write report output to"));

    CommandLine cmdline = new DefaultParser().parse(options, args);
    File annotationsDir = new File(cmdline.getOptionValue("annotationsDir"));
    File callFile       = new File(cmdline.getOptionValue("callFile"));
    Path outputFile     = Paths.get(cmdline.getOptionValue("reportFile"));

    new Reporter(annotationsDir)
        .analyze(callFile)
        .printMarkdown(outputFile);
  }

  /**
   * public constructor. start a new reporter based on annotation data found in the given <code>annotationsDir</code>.
   * @param annotationsDir directory of annotations files
   */
  public Reporter(@Nonnull File annotationsDir)  throws IOException {
    Preconditions.checkNotNull(annotationsDir);
    Preconditions.checkArgument(annotationsDir.exists());
    Preconditions.checkArgument(annotationsDir.isDirectory());

    File[] annotationFiles = annotationsDir.listFiles();
    if (annotationFiles == null || annotationFiles.length == 0) {
      throw new IOException("No annotation definitions to read from");
    }

    m_annotationFiles = Arrays.stream(annotationFiles)
        .filter(f -> f.getName().endsWith(".json"))
        .collect(Collectors.toList());
  }

  /**
   * Run the actual report process. Parse the input file, do the matching, and write the report files.
   * @param callFile file of haplotype calls
   */
  public Reporter analyze(@Nonnull File callFile) throws Exception {
    Preconditions.checkNotNull(callFile);
    Preconditions.checkArgument(callFile.exists());
    Preconditions.checkArgument(callFile.isFile());

    //Generate class used for loading JSON into
    JsonFileLoader loader = new JsonFileLoader();

    //Load the haplotype json, this is pointed at a test json and will likely break when meeting real
    // requiring some if not all rewriting
    List<GeneCall> calls = loader.loadHaplotypeGeneCalls(callFile.toPath());

    //Load the gene drug interaction list. This currently only handles single gene drug m_guidelineFiles and will require updating to handle multi gene drug interaction
    List<DosingGuideline> guidelines = loader.loadGuidelines(m_annotationFiles);

    //This is the primary work flow for generating the report where calls are matched to exceptions and drug gene m_guidelineFiles based on reported haplotypes
    m_dataUnifier = new DataUnifier(calls, guidelines);

    return this;
  }

  /**
   * Print a Markdown file of compiled report data
   * @param reportFile directory to write output to
   */
  public void printMarkdown(@Nonnull Path reportFile) throws IOException {
    Preconditions.checkNotNull(reportFile);
    sf_logger.debug("Writing output to {}", reportFile);

    new MarkdownWriter(reportFile)
        .print(m_dataUnifier);
  }

  @Nullable
  public List<GuidelineReport> getGuidelineReports() {
    if (m_dataUnifier == null) {
      return null;
    }
    return m_dataUnifier.getGuidelineResults();
  }
}

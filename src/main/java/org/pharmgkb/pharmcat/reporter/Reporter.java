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
import com.google.common.base.Preconditions;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.pharmgkb.pharmcat.haplotype.model.GeneCall;
import org.pharmgkb.pharmcat.reporter.io.JsonFileLoader;
import org.pharmgkb.pharmcat.reporter.io.ReporterWriter;
import org.pharmgkb.pharmcat.reporter.model.DosingGuideline;
import org.pharmgkb.pharmcat.reporter.resultsJSON.Interaction;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


/**
 * This is contains the main function for running the CPIC reporting tool.
 *
 * As of today  (4-23-2016) this tool is still in development currently this tool 
 * there are many broken, missing, and commented out pieces of this code due to
 * ongoing work to produce a working output.  
 *
 * Please contact me (Greyson Twist) on slack or at gtwist@cmh.edu if there are 
 * questions and I will try to improve documentation to make thing more clear
 *
 *
 * @author greytwist
 *
 */
public class Reporter {
  private static final Logger sf_logger = LoggerFactory.getLogger(MethodHandles.lookup().lookupClass());

  private List<File> m_annotationFiles = null;
  private File m_callFile = null;
  private Path m_reportDir = null;


  /**
   * main
   * @param args command line args
   * @throws Exception
   */
  public static void main(String[] args) throws Exception {
    Options options = new Options();

    options.addOption(new Option("annotationsDir", true, "required - directory holding all the annotations files"));
    options.addOption(new Option("callFile", true, "required - file from the Haplotyper"));
    options.addOption(new Option("reportDir", true, "required - directory to store output reports"));

    CommandLine cmdline = new DefaultParser().parse(options, args);
    File annotationsDir = new File(cmdline.getOptionValue("annotationsDir"));
    File callFile       = new File(cmdline.getOptionValue("callFile"));
    Path outputDir      = Paths.get(cmdline.getOptionValue("reportDir"));

    //if minimal required parameters are set parse command line
    Reporter report = new Reporter(annotationsDir, callFile, outputDir);
    //run reporter workflow
    report.run();
  }

  /**
   * parse command line options
   * @param annotationsDir directory of annotations files
   * @param callFile file of haplotype calls
   * @param reportDir directory to write output to
   * @throws IOException
   */
  public Reporter(@Nonnull File annotationsDir, @Nonnull File callFile, @Nonnull Path reportDir)  throws IOException {
    Preconditions.checkNotNull(annotationsDir);
    Preconditions.checkArgument(annotationsDir.exists());
    Preconditions.checkArgument(annotationsDir.isDirectory());

    Preconditions.checkNotNull(reportDir);
    Preconditions.checkArgument(reportDir.toFile().exists());
    Preconditions.checkArgument(reportDir.toFile().isDirectory());

    Preconditions.checkNotNull(callFile);
    Preconditions.checkArgument(callFile.exists());
    Preconditions.checkArgument(callFile.isFile());

    File[] annotationFiles = annotationsDir.listFiles();
    if (annotationFiles == null || annotationFiles.length == 0) {
      throw new IOException("No annotation definitions to read from");
    }

    m_annotationFiles = Arrays.stream(annotationFiles)
        .filter(f -> f.getName().endsWith(".json"))
        .collect(Collectors.toList());

    m_reportDir = reportDir;
    sf_logger.debug("Writing output to {}", m_reportDir);

    m_callFile = callFile;
  }

  public void run() throws Exception {

    //Generate class used for loading JSON into
    JsonFileLoader loader = new JsonFileLoader();

    //Load the haplotype json, this is pointed at a test json and will likely break when meeting real
    // requiring some if not all rewriting
    List<GeneCall> calls = loader.loadHaplotypeGeneCalls(m_callFile.toPath());

    //Load the gene drug interaction list. This currently only handles single gene drug m_guidelineFiles and will require updating to handle multi gene drug interaction
    List<DosingGuideline> guidelines = loader.loadGuidelines(m_annotationFiles);

    //This is the primary work flow for generating the report where calls are matched to exceptions and drug gene m_guidelineFiles based on reported haplotypes
    DataUnifier checker = new DataUnifier(calls, guidelines); // prime with data
    List<Interaction> results = checker.findMatches(); // run the actual comparison

    new ReporterWriter(m_reportDir)
        .print(results, checker.getSymbolToGeneReportMap());

    sf_logger.info("Complete");
  }
}

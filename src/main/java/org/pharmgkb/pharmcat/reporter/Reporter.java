package org.pharmgkb.pharmcat.reporter;

import java.io.BufferedWriter;
import java.io.IOException;
import java.lang.invoke.MethodHandles;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import javax.annotation.Nonnull;
import javax.annotation.Nullable;
import com.github.jknack.handlebars.Handlebars;
import com.github.jknack.handlebars.Template;
import com.github.jknack.handlebars.helper.StringHelpers;
import com.github.jknack.handlebars.io.ClassPathTemplateLoader;
import com.google.common.base.Preconditions;
import org.pharmgkb.common.io.util.CliHelper;
import org.pharmgkb.pharmcat.annotation.AnnotationReader;
import org.pharmgkb.pharmcat.haplotype.model.GeneCall;
import org.pharmgkb.pharmcat.reporter.handlebars.ReportHelpers;
import org.pharmgkb.pharmcat.reporter.io.AstrolabeOutputParser;
import org.pharmgkb.pharmcat.reporter.io.JsonFileLoader;
import org.pharmgkb.pharmcat.reporter.io.ReportData;
import org.pharmgkb.pharmcat.reporter.model.AstrolabeCall;
import org.pharmgkb.pharmcat.reporter.model.GuidelinePackage;
import org.pharmgkb.pharmcat.reporter.model.PharmcatException;
import org.pharmgkb.pharmcat.reporter.model.result.GuidelineReport;


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
  private static final String sf_templateName = "report";
  private static final String sf_templatePrefix = "/org/pharmgkb/pharmcat/reporter";
  private List<Path> m_annotationFiles = null;
  private List<PharmcatException> m_exceptions = null;
  private ReportContext m_reportContext = null;

  /**
   * main
   * @param args command line args
   */
  public static void main(String[] args) {

    CliHelper cliHelper = new CliHelper(MethodHandles.lookup().lookupClass())
        .addOption("d", "annotations-dir", "directory of allele definition files", true, "d")
        .addOption("c", "call-file", "named allele call JSON file", true, "c")
        .addOption("e", "exceptions-file", "exceptions logic file", true, "e")
        .addOption("a", "astrolabe-file", "optional, astrolabe call file", false, "a")
        .addOption("o", "output-file", "file path to write HTML report to", true, "o")
        ;

    try {
      if (!cliHelper.parse(args)) {
        System.exit(1);
      }

      Path annotationsDir = cliHelper.getValidDirectory("d", false);
      Path callFile = cliHelper.getValidFile("c", true);
      Path exceptionsFile = cliHelper.getValidFile("e", true);
      Path astrolabeFile = cliHelper.hasOption("a") ? cliHelper.getValidFile("a", true) : null;
      Path outputFile = cliHelper.getPath("o");

      new Reporter(annotationsDir, exceptionsFile)
          .analyze(callFile, astrolabeFile)
          .printHtml(outputFile);

    } catch (Exception ex) {
      ex.printStackTrace();
    }
  }

  /**
   * public constructor. start a new reporter based on annotation data found in the given <code>annotationsDir</code>.
   *
   * @param annotationsDir directory of annotation files
   */
  public Reporter(@Nonnull Path annotationsDir, @Nonnull Path exceptionsFile) throws IOException {

    Preconditions.checkNotNull(annotationsDir);
    Preconditions.checkArgument(Files.exists(annotationsDir));
    Preconditions.checkArgument(Files.isDirectory(annotationsDir));

    Preconditions.checkNotNull(exceptionsFile);
    Preconditions.checkArgument(Files.exists(exceptionsFile));
    Preconditions.checkArgument(Files.isRegularFile(exceptionsFile));

    AnnotationReader annotationReader = new AnnotationReader();
    annotationReader.read(annotationsDir);

    m_annotationFiles = Files.list(annotationsDir)
        .filter(f -> f.getFileName().toString().endsWith(".json"))
        .collect(Collectors.toList());
    if (m_annotationFiles.size() == 0) {
      throw new IOException("No annotation definitions to read from");
    }

    m_exceptions = Files.lines(exceptionsFile)
        .skip(1) // skip the header
        .map(PharmcatException::new)
        .collect(Collectors.toList());

  }

  /**
   * Run the actual report process. Parse the input file, do the matching, and write the report files.
   *
   * @param callFile file of haplotype calls
   */
  public Reporter analyze(@Nonnull Path callFile, @Nullable Path astrolabeFile) throws Exception {
    Preconditions.checkNotNull(callFile);
    Preconditions.checkArgument(Files.exists(callFile));
    Preconditions.checkArgument(Files.isRegularFile(callFile));

    //Generate class used for loading JSON into
    JsonFileLoader loader = new JsonFileLoader();

    //Load the haplotype json, this is pointed at a test json and will likely break when meeting real
    // requiring some if not all rewriting
    List<GeneCall> calls = loader.loadHaplotypeGeneCalls(callFile);

    //Load the astrolabe calls if it's available
    List<AstrolabeCall> astrolabeCalls = new ArrayList<>();
    if (astrolabeFile != null) {
      astrolabeCalls = AstrolabeOutputParser.parse(astrolabeFile);
    }

    //Load the gene drug interaction list. This currently only handles single gene drug m_guidelineFiles and will require updating to handle multi gene drug interaction
    List<GuidelinePackage> guidelines = loader.loadGuidelines(m_annotationFiles);

    //This is the primary work flow for generating the report where calls are matched to exceptions and drug gene m_guidelineFiles based on reported haplotypes
    m_reportContext = new ReportContext(calls, astrolabeCalls, guidelines, m_exceptions);

    return this;
  }

  /**
   * Print a HTML file of compiled report data
   * @param reportFile file to write output to
   */
  protected void printHtml(@Nonnull Path reportFile) throws IOException {

    Map<String,Object> reportData = ReportData.compile(m_reportContext);

    Handlebars handlebars = new Handlebars(new ClassPathTemplateLoader(sf_templatePrefix));
    StringHelpers.register(handlebars);
    handlebars.registerHelpers(ReportHelpers.class);
    Template template = handlebars.compile(sf_templateName);
    String html = template.apply(reportData);

    try (BufferedWriter writer = Files.newBufferedWriter(reportFile, StandardCharsets.UTF_8)) {
      writer.write(html);
    }
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

  protected ReportContext getContext() {
    return m_reportContext;
  }
}

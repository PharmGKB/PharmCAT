package org.pharmgkb.pharmcat.reporter;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.lang.invoke.MethodHandles;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import javax.annotation.Nonnull;
import javax.annotation.Nullable;
import com.github.jknack.handlebars.Handlebars;
import com.github.jknack.handlebars.helper.StringHelpers;
import com.github.jknack.handlebars.io.ClassPathTemplateLoader;
import com.google.common.base.Preconditions;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import org.pharmgkb.common.io.util.CliHelper;
import org.pharmgkb.common.util.PathUtils;
import org.pharmgkb.pharmcat.haplotype.model.GeneCall;
import org.pharmgkb.pharmcat.reporter.handlebars.ReportHelpers;
import org.pharmgkb.pharmcat.reporter.io.JsonFileLoader;
import org.pharmgkb.pharmcat.reporter.io.OutsideCallParser;
import org.pharmgkb.pharmcat.reporter.model.GuidelinePackage;
import org.pharmgkb.pharmcat.reporter.model.MessageAnnotation;
import org.pharmgkb.pharmcat.reporter.model.OutsideCall;
import org.pharmgkb.pharmcat.util.DataManager;


/**
 * This is the main class for running the reporting tool. It's responsible for taking input of all the
 * necessary data files, parsing them, and running the reporter components.
 *
 * This can be run both on the command line and procedurally.
 *
 * @author greytwist
 * @author Ryan Whaley
 */
public class Reporter {
  private static final String FINAL_REPORT      = "report";
  private static final String sf_templatePrefix = "/org/pharmgkb/pharmcat/reporter";
  private static final String sf_messagesFile   = "org/pharmgkb/pharmcat/reporter/messages.json";

  private static final Gson sf_gson = new GsonBuilder().serializeNulls().excludeFieldsWithoutExposeAnnotation()
      .setPrettyPrinting().create();
  private List<Path> m_annotationFiles;
  private List<MessageAnnotation> m_messages;
  private ReportContext m_reportContext = null;

  /**
   * Main CLI
   * @param args command line args
   */
  public static void main(String[] args) {

    CliHelper cliHelper = new CliHelper(MethodHandles.lookup().lookupClass())
        .addOption("c", "call-file", "named allele call JSON file", true, "c")
        .addOption("a", "outside-call-file", "optional, outside call TSV file", false, "a")
        .addOption("o", "output-file", "file path to write HTML report to", true, "o")
        .addOption("t", "title", "optional, text to add to the report title", false, "t")
        .addOption("g", "guidelines-dir", "directory of guideline annotations (JSON files)", false, "n")
        ;

    try {
      if (!cliHelper.parse(args)) {
        System.exit(1);
      }

      Path guidelinesDir = null;
      if (cliHelper.hasOption("g")) {
        guidelinesDir = cliHelper.getValidDirectory("g", false);
      }
      Path callFile = cliHelper.getValidFile("c", true);
      Path outsideCallPath = cliHelper.hasOption("a") ? cliHelper.getValidFile("a", true) : null;
      Path outputFile = cliHelper.getPath("o");
      String title = cliHelper.getValue("t");

      new Reporter(guidelinesDir)
          .analyze(callFile, outsideCallPath)
          .printHtml(outputFile, title, null);

    } catch (Exception ex) {
      ex.printStackTrace();
    }
  }

  /**
   * public constructor. start a new reporter based on annotation data found in the given <code>annotationsDir</code>.
   *
   * @param annotationsDir directory of annotation files
   */
  public Reporter(@Nullable Path annotationsDir) throws IOException {

    if (annotationsDir == null) {
      annotationsDir = DataManager.DEFAULT_GUIDELINE_DIR;
    }
    Preconditions.checkArgument(Files.exists(annotationsDir));
    Preconditions.checkArgument(Files.isDirectory(annotationsDir));

    m_annotationFiles = Files.list(annotationsDir)
        .filter(f -> f.getFileName().toString().endsWith(".json"))
        .collect(Collectors.toList());
    if (m_annotationFiles.size() == 0) {
      throw new IOException("No annotation definitions to read from");
    }

    try (BufferedReader reader = Files.newBufferedReader(PathUtils.getPathToResource(sf_messagesFile))) {
      MessageAnnotation[] messages = new Gson().fromJson(reader, MessageAnnotation[].class);
      m_messages = Arrays.asList(messages);
    }
  }

  /**
   * Run the actual report process. Parse the input file, do the matching, and write the report files.
   *
   * @param callFile file of haplotype calls
   */
  public Reporter analyze(@Nullable Path callFile, @Nullable Path outsideCallPath) throws Exception {
    //Generate class used for loading JSON into
    JsonFileLoader loader = new JsonFileLoader();

    //Load the haplotype json, this is pointed at a test json and will likely break when meeting real
    // requiring some if not all rewriting
    List<GeneCall> calls = new ArrayList<>();
    if (callFile != null) {
      Preconditions.checkArgument(Files.exists(callFile));
      Preconditions.checkArgument(Files.isRegularFile(callFile));
      calls = loader.loadHaplotypeGeneCalls(callFile);
    }

    //Load the outside calls if it's available
    List<OutsideCall> outsideCalls = new ArrayList<>();
    if (outsideCallPath != null) {
      Preconditions.checkArgument(Files.exists(outsideCallPath));
      Preconditions.checkArgument(Files.isRegularFile(outsideCallPath));
      outsideCalls = OutsideCallParser.parse(outsideCallPath);
    }

    //Load the gene drug interaction list. This currently only handles single gene-drug guidelines and will require updating to handle multi gene-drug interaction
    List<GuidelinePackage> guidelines = loader.loadGuidelines(m_annotationFiles);

    //This is the primary work flow for generating the report where calls are matched to exceptions and drug gene m_guidelineFiles based on reported haplotypes
    m_reportContext = new ReportContext(calls, outsideCalls, guidelines);

    m_reportContext.applyMessage(m_messages);

    return this;
  }

  /**
   * Print a HTML file of compiled report data
   * @param reportFile file to write output to
   */
  public void printHtml(@Nonnull Path reportFile, @Nullable String title, @Nullable Path jsonFile) throws IOException {

    Map<String,Object> reportData = m_reportContext.compile(title);

    Handlebars handlebars = new Handlebars(new ClassPathTemplateLoader(sf_templatePrefix));
    StringHelpers.register(handlebars);
    handlebars.registerHelpers(ReportHelpers.class);

    try (BufferedWriter writer = Files.newBufferedWriter(reportFile, StandardCharsets.UTF_8)) {
      writer.write(handlebars.compile(FINAL_REPORT).apply(reportData));
    }

    if (jsonFile != null) {
      try (BufferedWriter writer = Files.newBufferedWriter(jsonFile, StandardCharsets.UTF_8)) {
        writer.write(sf_gson.toJson(reportData));
        System.out.println("Writing JSON to " + jsonFile);
      }
    }
  }

  public ReportContext getContext() {
    return m_reportContext;
  }
}

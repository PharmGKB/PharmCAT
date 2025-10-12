package org.pharmgkb.pharmcat.stats;

import java.io.IOException;
import java.lang.invoke.MethodHandles;
import java.nio.file.DirectoryStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import com.google.common.base.Stopwatch;
import org.jspecify.annotations.Nullable;
import org.pharmgkb.common.util.CliHelper;
import org.pharmgkb.common.util.TimeUtils;
import org.pharmgkb.pharmcat.BaseConfig;
import org.pharmgkb.pharmcat.Env;
import org.pharmgkb.pharmcat.ReportableException;
import org.pharmgkb.pharmcat.UnexpectedStateException;
import org.pharmgkb.pharmcat.haplotype.ResultSerializer;
import org.pharmgkb.pharmcat.haplotype.model.Metadata;
import org.pharmgkb.pharmcat.phenotype.Phenotyper;
import org.pharmgkb.pharmcat.reporter.ReportContext;
import org.pharmgkb.pharmcat.reporter.format.CallsOnlyFormat;
import org.pharmgkb.pharmcat.util.CliUtils;


/**
 * This utility generates and merges calls-only TSV based on *.match.json files.
 *
 * @author Mark Woon
 */
public class MergeCalls {
  private static final Pattern sf_fileIdPattern = Pattern.compile("^(?:.+\\.)?(\\d+)\\.match\\.json");
  private enum MergeMode { MERGED_FILE, MERGED_FILE_PER_SUBDIRECTORY, INDIVIDUAL_FILES }
  private final Env m_env;
  private final @Nullable Path m_sampleMetadataFile;
  private final NumberFormat m_numFormat = NumberFormat.getInstance();


  public static void main(String[] args) {
    try {
      CliHelper cliHelper = new CliHelper(MethodHandles.lookup().lookupClass())
          .addOption("i", "input-dir", "Input directory", true, "directory")
          .addOption("o", "output-dir", "Output to (optional, default is input file directory)", false, "directory")
          .addOption("o1d", "output-1-per-subdir", "Output 1 file per subdirectory")
          .addOption("o1f", "output-1-per-file", "Output 1 file per *.match.json file")
          .addOption("sm", "sample-metadata", "TSV file containing sample metadata", false, "file")
          .addOption("esid", "extract-sample-id-from-filename", "Attempt to extract Sample ID from filename")
          .addOption("debug", "debug", "Add debug information to output files")
          ;

      if (!cliHelper.parse(args)) {
        if (!cliHelper.isHelpRequested() && !cliHelper.isVersionRequested()) {
          CliUtils.failIfNotTest();
        }
        return;
      }
      if (cliHelper.hasOption("o1d") && cliHelper.hasOption("o1f")) {
        System.out.println("-o1d and -o1f are mutually exclusive! Exiting...");
        System.exit(1);
      }


      Path inDir = cliHelper.getValidDirectory("i", false);
      Path outDir;
      if (cliHelper.hasOption("o")) {
        outDir = cliHelper.getValidDirectory("o", true);
      } else {
        outDir = inDir;
      }
      Path sampleMetadataFile = null;
      if (cliHelper.hasOption("sm")) {
        sampleMetadataFile = cliHelper.getValidFile("sm", true);
      }

      MergeMode mergeMode = MergeMode.MERGED_FILE;
      if (cliHelper.hasOption("o1d")) {
        mergeMode = MergeMode.MERGED_FILE_PER_SUBDIRECTORY;
      } else if (cliHelper.hasOption("o1f")) {
        mergeMode = MergeMode.INDIVIDUAL_FILES;
      }

      if (cliHelper.hasOption("debug")) {
        System.setProperty(CallsOnlyFormat.ENV_DEBUG_KEY, "true");
      }

      MergeCalls summarizer = new MergeCalls(sampleMetadataFile);
      summarizer.mergeFromMatchFiles(inDir, outDir, mergeMode, cliHelper.hasOption("esid"));

    } catch (CliHelper.InvalidPathException ex) {
      System.out.println(ex.getMessage());
    } catch (Exception e) {
      //noinspection CallToPrintStackTrace
      e.printStackTrace();
    }
  }


  private MergeCalls(@Nullable Path sampleMetadataFile) throws IOException, ReportableException {
    m_env = new Env();
    m_sampleMetadataFile = sampleMetadataFile;
  }


  /**
   * Generates merged calls-only TSV file(s) based on *.match.json files.
   *
   * @param inDir - directory to search recursively for *.match.json files
   * @param outDir - directory to save results to
   * @param lookupSampleIdFromFilename - true to attempt to get Sample ID from filename if not present
   */
  private void mergeFromMatchFiles(Path inDir, Path outDir, MergeMode mergeMode, boolean lookupSampleIdFromFilename)
      throws IOException {
    System.out.println("Reading from " + inDir);
    System.out.println("Writing to " + outDir);

    if (!Files.exists(outDir)) {
      Files.createDirectories(outDir);
    }

    readDir(inDir, outDir,
        mergeMode, lookupSampleIdFromFilename);
  }


  private void readDir(Path inDir, Path outDir, MergeMode mergeMode, boolean lookupSampleIdFromFilename)
      throws IOException {

    // CODE FOR SKIPPING NUMBERED DIRECTORIES:
//    try {
//      int subdir = Integer.parseInt(inDir.getFileName().toString());
//      if (subdir < 34) {
//        return;
//      }
//    } catch (NumberFormatException ex) {
//      // ignore
//    }

    System.out.println();
    System.out.println("Looking in "+ inDir);
    List<Path> subDirs = new ArrayList<>();

    int numFiles = 0;
    Stopwatch stopwatch = Stopwatch.createStarted();
    try (DirectoryStream<Path> stream = Files.newDirectoryStream(inDir)) {
      for (Path path : stream) {
        if (Files.isDirectory(path)) {
          subDirs.add(path);
        } else if (Files.isRegularFile(path) && path.toString().endsWith(".match.json")) {
          numFiles += 1;
          if (numFiles % 1000 == 0) {
            System.out.println(m_numFormat.format(numFiles) + " - " + TimeUtils.humanReadableDuration(stopwatch.elapsed()));
            stopwatch.reset()
                .start();
          }

          handleFile(inDir, outDir, mergeMode, lookupSampleIdFromFilename, path);
        }
      }
    }
    if (numFiles > 0 && numFiles % 1000 != 0) {
      System.out.println(m_numFormat.format(numFiles) + " - " + TimeUtils.humanReadableDuration(stopwatch.elapsed()));
    }
    System.out.println("---");

    for (Path d : subDirs) {
      readDir(d, outDir, mergeMode, lookupSampleIdFromFilename);
    }
  }

  private void handleFile(Path inDir, Path outDir, MergeMode mergeMode, boolean lookupSampleIdFromFilename,
      Path matchFile) throws IOException {

    org.pharmgkb.pharmcat.haplotype.model.Result matcherResult = new ResultSerializer()
        .fromJson(matchFile);

    Metadata metadata = matcherResult.getMetadata();
    // sample id was not available in results generated before 3.0
    //noinspection ConstantValue
    if (lookupSampleIdFromFilename && metadata.getSampleId() == null) {
      // add sample id based on filename pattern
      Matcher m = sf_fileIdPattern.matcher(matchFile.getFileName().toString());
      if (m.matches()) {
        metadata.setSampleId(m.group(1));
      }
    }
    //noinspection ConstantValue
    if (m_sampleMetadataFile != null && metadata.getSampleId() != null) {
      // add sample metadata
      Map<String, String> sampleData = m_env.getSampleMetadata(m_sampleMetadataFile, metadata.getSampleId(), true);
      if (sampleData != null && !sampleData.isEmpty()) {
        metadata.setSampleProps(sampleData);
      }
    }

    // don't rely on any existing phenotype.json file because it may not have the data we expect
    Phenotyper phenotyper = new Phenotyper(m_env, matcherResult.getMetadata(), matcherResult.getGeneCalls(),
        new TreeSet<>(), new HashMap<>());
    try {
      Path outFile;
      if (mergeMode == MergeMode.INDIVIDUAL_FILES) {
        // write tsv report next to matcher.json
        outFile = inDir.resolve(BaseConfig.getBaseFilename(matchFile) + BaseConfig.REPORTER_SUFFIX + ".tsv");
      } else {
        // write tsv report to file in the output directory
        if (mergeMode == MergeMode.MERGED_FILE) {
          outFile = outDir.resolve("merged" + BaseConfig.REPORTER_SUFFIX + ".tsv");
        } else {
          outFile = outDir.resolve("merged_" + inDir.getFileName() + BaseConfig.REPORTER_SUFFIX + ".tsv");
        }
      }

      ReportContext reportContext = new ReportContext(m_env, phenotyper, BaseConfig.getBaseFilename(matchFile));
      CallsOnlyFormat callsOnlyFormat = new CallsOnlyFormat(outFile, m_env);
      if (mergeMode != MergeMode.INDIVIDUAL_FILES) {
        callsOnlyFormat.singleFileMode();
      }
      callsOnlyFormat.write(reportContext);

    } catch (UnexpectedStateException ex) {
      System.out.println("Error with " + matchFile);
      //noinspection CallToPrintStackTrace
      ex.printStackTrace();
    }
  }
}

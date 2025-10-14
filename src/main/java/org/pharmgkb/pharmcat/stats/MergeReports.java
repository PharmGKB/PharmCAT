package org.pharmgkb.pharmcat.stats;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.lang.invoke.MethodHandles;
import java.nio.file.DirectoryStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import com.google.common.base.Splitter;
import com.google.common.base.Stopwatch;
import org.apache.commons.compress.archivers.tar.TarArchiveEntry;
import org.apache.commons.compress.archivers.tar.TarArchiveInputStream;
import org.apache.commons.compress.compressors.gzip.GzipCompressorInputStream;
import org.apache.commons.lang3.StringUtils;
import org.jspecify.annotations.Nullable;
import org.pharmgkb.common.util.CliHelper;
import org.pharmgkb.common.util.StreamUtils;
import org.pharmgkb.common.util.TimeUtils;

import static org.pharmgkb.pharmcat.reporter.format.CallsOnlyFormat.HEADER_SAMPLE_ID;


/**
 * This utility merges report TSV's and adds metadata.
 *
 * @author Mark Woon
 */
public class MergeReports {
  // expect sample ID to be composed of letters and numbers
  private static final Pattern sf_fileIdPattern = Pattern.compile("^(?:.+\\.)?([A-Za-z\\d]+)\\.report\\.tsv");
  private static final Splitter sf_tabSplitter = Splitter.on("\t").trimResults();
  private final boolean m_anonymize;
  private final boolean m_verbose;
  private final NumberFormat m_numFormat = NumberFormat.getInstance();
  private final boolean m_lookupSampleIdFromFilename;
  private final @Nullable List<String> m_metadataHeaders;
  private final @Nullable Map<String, List<String>> m_metadata;
  private @Nullable Stopwatch m_stopwatch;
  private @Nullable String m_version;
  private int m_numFiles;
  private final List<String> m_missingMetadata = new ArrayList<>();


  public static void main(String[] args) {
    try {
      CliHelper cliHelper = new CliHelper(MethodHandles.lookup().lookupClass())
          .addOption("i", "input-dir", "Input directory", true, "directory")
          .addOption("o", "output-file",
              "Output file (optional, defaults to 'merged.reports.tsv' in input file directory)", false, "file")
          .addOption("sm", "sample-metadata", "TSV file containing sample metadata", false, "file")
          .addOption("esid", "extract-sample-id-from-filename", "Attempt to extract Sample ID from filename")
          .addOption("anon", "anonymize", "Hide Sample IDs")
          ;

      if (!cliHelper.parse(args)) {
        return;
      }

      Path inDir = cliHelper.getValidDirectory("i", false);
      Path outFile;
      if (cliHelper.hasOption("o")) {
        outFile = cliHelper.getValidFile("o", false);
      } else {
        outFile = inDir.resolve("merged.reports.tsv");
      }
      Path sampleMetadataFile = null;
      if (cliHelper.hasOption("sm")) {
        sampleMetadataFile = cliHelper.getValidFile("sm", true);
      }

      MergeReports summarizer = new MergeReports(sampleMetadataFile, cliHelper.hasOption("esid"),
          cliHelper.hasOption("anon"), cliHelper.isVerbose());
      summarizer.mergeFromReports(inDir, outFile);

    } catch (CliHelper.InvalidPathException ex) {
      System.out.println(ex.getMessage());
    } catch (Exception e) {
      //noinspection CallToPrintStackTrace
      e.printStackTrace();
    }
  }


  private MergeReports(@Nullable Path sampleMetadataFile, boolean lookupSampleIdFromFilename, boolean anonymize,
      boolean verbose) throws IOException {
    m_anonymize = anonymize;
    m_verbose = verbose;
    m_lookupSampleIdFromFilename = lookupSampleIdFromFilename;
    if (sampleMetadataFile != null) {
      m_metadata = new HashMap<>();
      try (BufferedReader reader = StreamUtils.openReader(sampleMetadataFile)) {
        // expect headers
        String line = reader.readLine();
        m_metadataHeaders = new ArrayList<>(sf_tabSplitter.splitToList(line));
        if (m_metadataHeaders.size() < 2) {
          throw new IllegalStateException("Expecting at least two columns in sample metadata file");
        }
        // delete the first header (Sample ID column)
        m_metadataHeaders.remove(0);

        while ((line = reader.readLine()) != null) {
          List<String> data = sf_tabSplitter.splitToList(line);
          m_metadata.put(data.get(0), data.subList(1, data.size()));
        }
      }
    } else {
      m_metadataHeaders = null;
      m_metadata = null;
    }
  }


  /**
   * Generates merged calls-only TSV file(s) based on *.report.tsv files.
   *
   * @param inDir - directory to search recursively for *.report.tsv files
   * @param outFile - file to save results to
   */
  private void mergeFromReports(Path inDir, Path outFile)
      throws IOException {
    System.out.println("Reading from " + inDir);
    System.out.println("Writing to " + outFile);

    if (m_verbose) {
      m_stopwatch = Stopwatch.createStarted();
    }
    try (PrintWriter writer = new PrintWriter(Files.newBufferedWriter(outFile))) {
      readDir(inDir, writer);
    }
  }

  private void incrementNumFiles() {
    m_numFiles += 1;
    if (m_numFiles % 1000 == 0) {
      printStatus();
    }
  }

  private void printStatus() {
    if (m_verbose && m_numFiles > 0) {
      System.out.print(m_numFormat.format(m_numFiles));
      if (m_stopwatch != null) {
        System.out.print(" - " + TimeUtils.humanReadableDuration(m_stopwatch.elapsed()));
        m_stopwatch.reset()
            .start();
      }
      System.out.println();
    }
  }


  private void readDir(Path inDir, PrintWriter writer)
      throws IOException {

    System.out.println();
    System.out.println("Looking in " + inDir);
    List<Path> subDirs = new ArrayList<>();

    try (DirectoryStream<Path> stream = Files.newDirectoryStream(inDir)) {
      for (Path path : stream) {
        if (Files.isDirectory(path)) {
          subDirs.add(path);
        } else if (Files.isRegularFile(path)) {
          if (path.toString().endsWith(".report.tsv")) {
            handleFile(path, writer);
          } else if ((path.toString().endsWith(".tgz") || path.toString().endsWith(".tar.gz"))) {
            readTgz(path, writer);
          }
        }
      }
    }
    printStatus();

    for (Path d : subDirs) {
      readDir(d, writer);
    }
  }

  private void handleFile(Path matchFile, PrintWriter writer) throws IOException {
    try (BufferedReader reader = Files.newBufferedReader(matchFile)) {
      handleFile(matchFile.getFileName().toString(), reader, writer);
    }
  }

  private void handleFile(String filename, BufferedReader reader, PrintWriter writer) throws IOException {
    // expect PharmCAT version
    String version = reader.readLine();
    if (!version.startsWith("PharmCAT")) {
      throw new IllegalStateException("Expecting PharmCAT version header but found: " + version);
    }
    if (m_version == null) {
      m_version = version;
    } else if (!m_version.equals(version)) {
      throw new IllegalStateException("Found different PharmCAT versions: " + m_version + " vs. " + version);
    }
    // expect headers
    String line = reader.readLine();
    List<String> origHeaders = sf_tabSplitter.splitToList(line);
    boolean fileHasSampleId = origHeaders.get(0).equals(HEADER_SAMPLE_ID);
    String sampleId = null;
    if (m_metadata != null) {
      // need Sample ID to look up metadata
      // an existing sample id column takes precedence
      if (m_lookupSampleIdFromFilename && !fileHasSampleId) {
        // add sample id based on filename pattern
        Matcher m = sf_fileIdPattern.matcher(filename);
        if (m.matches()) {
          sampleId = m.group(1);
        }
      }
    }

    List<String> newHeaders = new ArrayList<>(origHeaders);
    if (m_anonymize) {
      if (fileHasSampleId) {
        newHeaders.remove(HEADER_SAMPLE_ID);
      }
    } else {
      newHeaders.add(0, HEADER_SAMPLE_ID);
    }
    if (m_metadataHeaders != null) {
      newHeaders.addAll(m_metadataHeaders);
    }
    if (m_numFiles == 0) {
      // write headers
      writer.println(m_version);
      writer.println(String.join("\t", newHeaders));
    }

    while ((line = reader.readLine()) != null) {
      if (StringUtils.isBlank(line)) {
        continue;
      }
      List<String> data = sf_tabSplitter.splitToList(line);
      if (sampleId == null && fileHasSampleId) {
        sampleId = data.get(0);
      }

      if (m_anonymize) {
        if (fileHasSampleId) {
          data.remove(0);
        }
      } else {
        if (sampleId != null) {
          writer.print(sampleId);
        } else {
          System.out.println("Warning: no sample id for file " + filename);
        }
        writer.print("\t");
      }
      writer.print(String.join("\t", data));
      if (m_metadata != null && sampleId != null) {
        if (m_metadata.containsKey(sampleId)) {
          writer.print("\t");
          writer.print(String.join("\t", m_metadata.get(sampleId)));
        } else {
          if (!m_missingMetadata.contains(sampleId)) {
            m_missingMetadata.add(sampleId);
            System.out.println("Warning: no metadata for sample " + sampleId + " in " + filename);
          }
        }
      }
      writer.println();
    }
    incrementNumFiles();
  }

  /**
   * Reads a .tgz or .tar.gz archive and processes any contained *.report.tsv files.
   * Returns the number of files processed.
   */
  private void readTgz(Path tgzFile, PrintWriter writer) throws IOException {
    try (BufferedInputStream bis = new BufferedInputStream(Files.newInputStream(tgzFile));
         GzipCompressorInputStream gzis = new GzipCompressorInputStream(bis);
         TarArchiveInputStream tais = new TarArchiveInputStream(gzis)) {

      System.out.println("Looking in " + tgzFile);
      TarArchiveEntry entry;
      while ((entry = tais.getNextEntry()) != null) {
        if (entry.isDirectory()) {
          continue;
        }
        String name = entry.getName();
        if (name != null) {
          int idx = name.lastIndexOf('/');
          if (idx >= 0) {
            name = name.substring(idx + 1);
          }
          if (name.endsWith(".report.tsv")) {
            handleFile(name, new BufferedReader(new InputStreamReader(tais)), writer);
          }
        }
      }
    }
  }
}

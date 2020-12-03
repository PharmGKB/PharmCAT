package org.pharmgkb.pharmcat.util;

import java.io.BufferedReader;
import java.io.IOException;
import java.lang.invoke.MethodHandles;
import java.net.URI;
import java.nio.file.DirectoryStream;
import java.nio.file.FileSystem;
import java.nio.file.FileSystems;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.util.HashMap;
import java.util.Map;
import java.util.Properties;
import java.util.Set;
import java.util.stream.Collectors;
import com.google.common.base.Preconditions;
import org.apache.commons.io.FileUtils;
import org.apache.commons.lang3.StringUtils;
import org.pharmgkb.common.io.util.CliHelper;
import org.pharmgkb.common.util.PathUtils;
import org.pharmgkb.pharmcat.definition.CuratedDefinitionParser;
import org.pharmgkb.pharmcat.definition.model.DefinitionFile;


/**
 * This class manages external resources (e.g. allele definition files, dosing guideline annotations)
 *
 * @author Mark Woon
 */
public class DataManager {
  public static final Path DEFAULT_DEFINITION_DIR = PathUtils.getPathToResource("org/pharmgkb/pharmcat/definition/alleles");
  private static final Path DEFAULT_REPORTER_DIR = PathUtils.getPathToResource("org/pharmgkb/pharmcat/reporter");
  public static final String EXEMPTIONS_JSON_FILE_NAME = "exemptions.json";
  private static final String MESSAGES_JSON_FILE_NAME = "messages.json";
  private static final String ALLELE_DEFINITION_ARCHIVE =  "allele.definitions.zip";
  private final String m_googleUser;
  private final String m_googleKey;
  private final DataSerializer m_dataSerializer = new DataSerializer();
  private final boolean m_verbose;


  private DataManager(Path propertyFile, boolean verbose) throws IOException {

    Properties properties = new Properties();
    try (BufferedReader reader = Files.newBufferedReader(propertyFile)) {
      properties.load(reader);
    }
    m_googleUser = StringUtils.stripToNull((String)properties.get("google.user"));
    Preconditions.checkState(m_googleUser != null, "Missing property: 'google.user");
    m_googleKey = StringUtils.stripToNull((String)properties.get("google.key"));
    Preconditions.checkState(m_googleKey != null, "Missing property: 'google.key");
    m_verbose = verbose;
  }


  public static void main(String[] args) {

    try {
      CliHelper cliHelper = new CliHelper(MethodHandles.lookup().lookupClass())
          .addOption("p", "properties-file", "PharmCAT properties file", false, "p")
          .addOption("dl", "download-dir", "directory to save downloaded files", false, "dl")
          .addOption("a", "alleles-dir", "directory to save generated allele definition files", false, "a")
          .addOption("m", "messages-dir", "directory to write messages to", false, "m")
          .addOption("g", "guidelines-dir", "directory to save guideline annotations to", false, "g")
          .addOption("sd", "skip-download", "skip downloading")
          .addOption("sa", "skip-alleles", "skip alleles")
          .addOption("sm", "skip-messages", "skip messages")
          .addOption("sg", "skip-guidelines", "skip guidelines");


      if (!cliHelper.parse(args)) {
        System.exit(1);
      }

      Path propsFile = CliUtils.getPropsFile(cliHelper, "p");
      Path downloadDir;
      if (cliHelper.hasOption("dl")) {
        downloadDir = cliHelper.getValidDirectory("dl", true);
      } else {
        downloadDir = Files.createTempDirectory("pharmcat");
        downloadDir.toFile().deleteOnExit();
      }

      try {
        Path allelesDir;
        if (cliHelper.hasOption("a")) {
          allelesDir = cliHelper.getValidDirectory("a", true);
        } else {
          allelesDir = DEFAULT_DEFINITION_DIR;
        }
        Path messageDir;
        if (cliHelper.hasOption("m")) {
          messageDir = cliHelper.getValidDirectory("m", true);
        } else {
          messageDir = DEFAULT_REPORTER_DIR;
        }

        Path exemptionsTsv = downloadDir.resolve("exemptions.tsv");
        Path exemptionsJson = allelesDir.resolve(EXEMPTIONS_JSON_FILE_NAME);
        Path messagesTsv = downloadDir.resolve("messages.tsv");
        Path messagesJson = messageDir.resolve(MESSAGES_JSON_FILE_NAME);

        DataManager manager = new DataManager(propsFile, cliHelper.isVerbose());
        if (!cliHelper.hasOption("sd")) {
          manager.download(downloadDir, exemptionsTsv, messagesTsv);
        }

        if (!cliHelper.hasOption("sa")) {
          manager.transformAlleleDefinitions(downloadDir, allelesDir);
          manager.transformExemptions(exemptionsTsv, exemptionsJson);
        }
        if (!cliHelper.hasOption("sm")) {
          manager.transformMessages(messagesTsv, messagesJson);
        }

      } finally {
        if (!cliHelper.hasOption("dl")) {
          FileUtils.deleteQuietly(downloadDir.toFile());
        }
      }

    } catch (Exception ex) {
      ex.printStackTrace();
    }
  }


  private void download(Path downloadDir, Path exemptionsFile, Path messagesFile) throws Exception {

    SheetsHelper sh = new SheetsHelper(m_googleUser, m_googleKey);
    sh.downloadAlleleDefinitions(downloadDir);
    sh.downloadAlleleExemptionsFile(exemptionsFile);
    sh.downloadMessagesFile(messagesFile);
  }


  /**
   * Does the work for stepping through the files and applying the format.
   */
  private void transformAlleleDefinitions(Path downloadDir, Path definitionsDir) throws Exception {

    System.out.println();
    System.out.println("Saving allele definitions in " + definitionsDir.toString());
    Set<String> currentFiles = Files.list(definitionsDir)
        .map(PathUtils::getFilename)
        .collect(Collectors.toSet());

    try (FileSystem zipFs = makeZipFilesystem(definitionsDir, ALLELE_DEFINITION_ARCHIVE)) {
      try (DirectoryStream<Path> files = Files.newDirectoryStream(downloadDir, f -> f.toString().endsWith("_translation.tsv"))) {
        for (Path file : files) {
          if (m_verbose) {
            System.out.println("Parsing " + file);
          }
          CuratedDefinitionParser parser = new CuratedDefinitionParser(file);

          DefinitionFile definitionFile = parser.parse();
          if (!parser.getWarnings().isEmpty()) {
            System.out.println("Warnings for " + file);
            parser.getWarnings()
                .forEach(System.out::println);
          }

          Path internalPath = zipFs.getPath(PathUtils.getFilename(file));
          Files.copy(file, internalPath, StandardCopyOption.REPLACE_EXISTING);

          Path jsonFile = definitionsDir.resolve(PathUtils.getBaseFilename(file) + ".json");
          m_dataSerializer.serializeToJson(definitionFile, jsonFile);
          if (m_verbose) {
            System.out.println("Wrote " + jsonFile);
          }
          currentFiles.remove(PathUtils.getFilename(jsonFile));
        }
      }
    }
    deleteObsoleteFiles(definitionsDir, currentFiles);
  }

  private void transformExemptions(Path tsvFile, Path jsonFile) throws IOException {

    System.out.println();
    System.out.println("Saving exemptions to " + jsonFile.toString());
    m_dataSerializer.serializeToJson(m_dataSerializer.deserializeExemptionsFromTsv(tsvFile), jsonFile);
  }

  private void transformMessages(Path tsvFile, Path jsonFile) throws IOException {

    System.out.println();
    System.out.println("Saving messages to " + jsonFile.toString());
    m_dataSerializer.serializeToJson(m_dataSerializer.deserializeMessagesFromTsv(tsvFile), jsonFile);
  }


  private void deleteObsoleteFiles(Path dir, Set<String> obsoleteFilenames) {

    for (String filename : obsoleteFilenames) {
      Path file = dir.resolve(filename);
      System.out.println("Deleting obsolete file: " + file);
      FileUtils.deleteQuietly(file.toFile());
    }
  }

  private FileSystem makeZipFilesystem(Path directory, String filename) throws Exception {

    Path zipFile = directory.resolve(filename);
    System.out.println();
    System.out.println("Saving allele definition tables to "+zipFile.toString());

    Map<String,String> env = new HashMap<>();
    env.put("create", String.valueOf(Files.notExists(zipFile)));
    URI fileUri = zipFile.toUri();
    URI zipUri = new URI("jar:" + fileUri.getScheme(), fileUri.getPath(), null);
    return FileSystems.newFileSystem(zipUri, env);
  }
}

package org.pharmgkb.pharmcat.util;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
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
import java.util.zip.ZipFile;
import javax.annotation.Nonnull;
import com.google.common.base.Preconditions;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.IOUtils;
import org.apache.commons.lang3.StringUtils;
import org.apache.http.HttpEntity;
import org.apache.http.client.methods.CloseableHttpResponse;
import org.apache.http.client.methods.HttpGet;
import org.apache.http.impl.client.CloseableHttpClient;
import org.apache.http.impl.client.HttpClients;
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
  public static final Path DEFAULT_REPORTER_DIR = PathUtils.getPathToResource("org/pharmgkb/pharmcat/reporter");
  public static final Path DEFAULT_GUIDELINE_DIR = PathUtils.getPathToResource("org/pharmgkb/pharmcat/reporter/guidelines");
  public static final String DOSING_GUIDELINE_URL = "https://api.pharmgkb.org/v1/download/file/data/dosingGuidelines.extended.json.zip?ref=pharmcat";
  public static final String EXEMPTIONS_JSON_FILE_NAME = "exemptions.json";
  public static final String MESSAGES_JSON_FILE_NAME = "messages.json";
  public static final String GUIDELINE_TIMESTAMP_FILE_NAME = "timestamp.txt";
  public static final String ALLELE_DEFINITION_ARCHIVE =  "allele.definitions.zip";
  private String m_googleUser;
  private String m_googleKey;
  private DataSerializer m_dataSerializer = new DataSerializer();
  private boolean m_verbose;



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
        Path guidelinesDir;
        if (cliHelper.hasOption("g")) {
          guidelinesDir = cliHelper.getValidDirectory("g", true);
        } else {
          guidelinesDir = DEFAULT_GUIDELINE_DIR;
          if (Files.exists(guidelinesDir)) {
            if (!Files.isDirectory(guidelinesDir)) {
              System.out.println(guidelinesDir + " is not a directory");
              return;
            }
          } else {
            Files.createDirectories(guidelinesDir);
          }
        }

        Path exemptionsTsv = downloadDir.resolve("exemptions.tsv");
        Path exemptionsJson = allelesDir.resolve(EXEMPTIONS_JSON_FILE_NAME);
        Path messagesTsv = downloadDir.resolve("messages.tsv");
        Path messagesJson = messageDir.resolve(MESSAGES_JSON_FILE_NAME);
        Path guidelinesZip = downloadDir.resolve("guidelines.zip");

        DataManager manager = new DataManager(propsFile, cliHelper.isVerbose());
        if (!cliHelper.hasOption("sd")) {
          manager.download(downloadDir, exemptionsTsv, messagesTsv);
          if (!cliHelper.hasOption("sg")) {
            manager.downloadGuidelines(guidelinesZip);
          }
        }

        if (!cliHelper.hasOption("sa")) {
          manager.transformAlleleDefinitions(downloadDir, allelesDir);
          manager.transformExemptions(exemptionsTsv, exemptionsJson);
        }
        if (!cliHelper.hasOption("sm")) {
          manager.transformMessages(messagesTsv, messagesJson);
        }
        if (!cliHelper.hasOption("sg")) {
          manager.transformGuidelines(guidelinesZip, guidelinesDir);
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


  private void download(@Nonnull Path downloadDir, @Nonnull Path exemptionsFile, @Nonnull Path messagesFile) throws Exception {

    SheetsHelper sh = new SheetsHelper(m_googleUser, m_googleKey);
    sh.downloadAlleleDefinitions(downloadDir);
    sh.downloadAlleleExemptionsFile(exemptionsFile);
    sh.downloadMessagesFile(messagesFile);
  }


  /**
   * Does the work for stepping through the files and applying the format.
   */
  private void transformAlleleDefinitions(@Nonnull Path downloadDir, @Nonnull Path definitionsDir) throws Exception {

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

  private void transformExemptions(@Nonnull Path tsvFile, @Nonnull Path jsonFile) throws IOException {

    System.out.println();
    System.out.println("Saving exemptions to " + jsonFile.toString());
    m_dataSerializer.serializeToJson(m_dataSerializer.deserializeExemptionsFromTsv(tsvFile), jsonFile);
  }

  private void transformMessages(@Nonnull Path tsvFile, @Nonnull Path jsonFile) throws IOException {

    System.out.println();
    System.out.println("Saving messages to " + jsonFile.toString());
    m_dataSerializer.serializeToJson(m_dataSerializer.deserializeMessagesFromTsv(tsvFile), jsonFile);
  }


  private void downloadGuidelines(@Nonnull Path guidelinesZip) throws Exception {

    try (CloseableHttpClient client = HttpClients.createDefault()) {
      try (CloseableHttpResponse response = client.execute(new HttpGet(DOSING_GUIDELINE_URL))) {
        HttpEntity entity = response.getEntity();
        if (entity != null) {
          try (InputStream in = entity.getContent();
               OutputStream out = Files.newOutputStream(guidelinesZip)) {
            IOUtils.copy(in, out);
          }
        }
      }
    }
  }

  private void transformGuidelines(@Nonnull Path zipFile, @Nonnull Path guidelinesDir) throws IOException {

    System.out.println();
    System.out.println("Saving guidelines to " + guidelinesDir.toString());

    Set<String> currentFiles = Files.list(guidelinesDir)
        .map(PathUtils::getFilename)
        .collect(Collectors.toSet());
    try (ZipFile zip = new ZipFile(zipFile.toFile())) {
      // extract timestamp
      zip.stream()
          .filter(ze -> ze.getName().toLowerCase().matches("created_.*\\.txt"))
          .forEachOrdered(ze -> {
            try (InputStream inputStream = zip.getInputStream(ze);
                 OutputStream out = Files.newOutputStream(guidelinesDir.resolve(GUIDELINE_TIMESTAMP_FILE_NAME))) {
              if (m_verbose) {
                System.out.println("Extracting " + ze.getName());
              }
              IOUtils.copy(inputStream, out);
              currentFiles.remove(GUIDELINE_TIMESTAMP_FILE_NAME);
            } catch (IOException ex) {
              throw new RuntimeException("Error extracting " + ze.getName(), ex);
            }
          });
      // exstract CPIC .json files
      zip.stream()
          .filter(ze -> ze.getName().toLowerCase().matches("cpic_.*\\.json"))
          .forEachOrdered(ze -> {
            Path file = guidelinesDir.resolve(ze.getName());
            try (InputStream inputStream = zip.getInputStream(ze);
                 OutputStream out = Files.newOutputStream(file)) {
              if (m_verbose) {
                System.out.println("Extracting " + ze.getName());
              }
              IOUtils.copy(inputStream, out);
              currentFiles.remove(PathUtils.getFilename(file));
            } catch (IOException ex) {
              throw new RuntimeException("Error extracting to " + file, ex);
            }
          });
      // remove obsolete files
      deleteObsoleteFiles(guidelinesDir, currentFiles);
    }
  }


  private void deleteObsoleteFiles(@Nonnull Path dir, @Nonnull Set<String> obsoleteFilenames) {

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

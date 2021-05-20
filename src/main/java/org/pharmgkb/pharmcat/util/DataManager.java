package org.pharmgkb.pharmcat.util;

import java.io.IOException;
import java.lang.invoke.MethodHandles;
import java.net.URL;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;
import com.google.common.base.Charsets;
import com.google.gson.Gson;
import org.apache.commons.io.FileUtils;
import org.pharmgkb.common.io.util.CliHelper;
import org.pharmgkb.common.util.PathUtils;
import org.pharmgkb.pharmcat.definition.model.DefinitionExemption;
import org.pharmgkb.pharmcat.definition.model.DefinitionFile;


/**
 * This class manages external resources (e.g. allele definition files, dosing guideline annotations)
 *
 * @author Mark Woon
 */
public class DataManager {
  public static final Path DEFAULT_DEFINITION_DIR = PathUtils.getPathToResource("org/pharmgkb/pharmcat/definition/alleles");
  public static final String EXEMPTIONS_JSON_FILE_NAME = "exemptions.json";
  private static final String MESSAGES_JSON_FILE_NAME = "messages.json";
  private final DataSerializer m_dataSerializer = new DataSerializer();
  private final boolean m_verbose;


  private DataManager(boolean verbose) {
    m_verbose = verbose;
  }


  public static void main(String[] args) {

    try {
      CliHelper cliHelper = new CliHelper(MethodHandles.lookup().lookupClass())
          .addOption("dl", "download-dir", "directory to save downloaded files", false, "dl")
          .addOption("a", "alleles-dir", "directory to save generated allele definition files", true, "a")
          .addOption("m", "messages-dir", "directory to write messages to", true, "m")
          .addOption("g", "guidelines-dir", "directory to save guideline annotations to", false, "g")
          .addOption("sd", "skip-download", "skip downloading")
          .addOption("sa", "skip-alleles", "skip alleles")
          .addOption("sm", "skip-messages", "skip messages")
          .addOption("sg", "skip-guidelines", "skip guidelines");


      if (!cliHelper.parse(args)) {
        System.exit(1);
      }

      Path downloadDir;
      if (cliHelper.hasOption("dl")) {
        downloadDir = cliHelper.getValidDirectory("dl", true);
      } else {
        downloadDir = Files.createTempDirectory("pharmcat");
        downloadDir.toFile().deleteOnExit();
      }

      try {
        Path allelesDir = cliHelper.getValidDirectory("a", true);
        Path messageDir = cliHelper.getValidDirectory("m", true);

        Path guidelineDir = null;
        if (cliHelper.hasOption("g")) {
          guidelineDir = cliHelper.getValidDirectory("g", true);
        }

        Path exemptionsTsv = downloadDir.resolve("exemptions.tsv");
        Path exemptionsJson = allelesDir.resolve(EXEMPTIONS_JSON_FILE_NAME);
        Path messagesTsv = downloadDir.resolve("messages.tsv");
        Path messagesJson = messageDir.resolve(MESSAGES_JSON_FILE_NAME);

        DataManager manager = new DataManager(cliHelper.isVerbose());
        if (!cliHelper.hasOption("sd")) {
          manager.download(downloadDir, exemptionsTsv, messagesTsv, guidelineDir);
        }

        if (!cliHelper.hasOption("sa")) {
          Map<String, DefinitionExemption> exemptionsMap = manager.transformExemptions(exemptionsTsv, exemptionsJson);
          manager.transformAlleleDefinitions(downloadDir, allelesDir, exemptionsMap);
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


  private void download(Path downloadDir, Path exemptionsFile, Path messagesFile, Path guidelinesDir) throws Exception {

    // download allele definitions
    FileUtils.copyURLToFile(
        new URL("https://files.cpicpgx.org/data/report/current/allele_definitions.json"),
        downloadDir.resolve("allele_definitions.json").toFile());

    // download drugs/guidelines file
    if (guidelinesDir != null) {
      Path drugsPath = guidelinesDir.resolve("drugs.json");
      FileUtils.copyURLToFile(
          new URL("https://files.cpicpgx.org/data/report/current/drugs.json"),
          drugsPath.toFile());
      System.out.println();
      System.out.println("Saving messages to " + drugsPath);    }

    String urlFmt = "https://docs.google.com/spreadsheets/d/%s/export?format=tsv";
    // download exemptions and messages
    FileUtils.copyURLToFile(new URL(String.format(urlFmt, "1xHvvXQIMv3xbqNhuN7zG6WP4DB7lpQDmLvz18w-u_lk")),
        exemptionsFile.toFile());
    FileUtils.copyURLToFile(new URL(String.format(urlFmt, "1MkWV6TlJTnw-KRNWeylyUJAUocCgupcJLlmV2fRdtcM")),
        messagesFile.toFile());
  }


  /**
   * Does the work for stepping through the files and applying the format.
   */
  private void transformAlleleDefinitions(Path downloadDir, Path definitionsDir,
      Map<String, DefinitionExemption> exemptionsMap) throws Exception {

    Path definitionsFile = downloadDir.resolve("allele_definitions.json");
    if (!Files.exists(definitionsFile)) {
      throw new IOException("Cannot find alleles definitions (" + definitionsFile + ")");
    }
    String json = FileUtils.readFileToString(definitionsFile.toFile(), Charsets.UTF_8);
    Gson gson = new Gson();
    Map<String, DefinitionFile> definitionFileMap = new HashMap<>();
    for (DefinitionFile df : gson.fromJson(json, DefinitionFile[].class)) {
      definitionFileMap.put(df.getGeneSymbol(), df);
    }

    System.out.println();
    System.out.println("Saving allele definitions in " + definitionsDir.toString());
    Set<String> currentFiles = Files.list(definitionsDir)
        .map(PathUtils::getFilename)
        .filter(f -> f.endsWith("_translation.json"))
        .collect(Collectors.toSet());

    for (String gene : exemptionsMap.keySet()) {
      DefinitionExemption exemption = exemptionsMap.get(gene);
      if (exemption.getIgnoredPositions().size() > 0) {
        System.out.println("Removing ignored positions in " + gene + "...");
        DefinitionFile definitionFile = definitionFileMap.get(gene);
        if (definitionFile == null) {
          throw new Exception("No definition for " + gene);
        }
        definitionFile.removeIgnoredPositions(exemption);
      }
    }

    for (String gene : definitionFileMap.keySet()) {
      DefinitionFile definitionFile = definitionFileMap.get(gene);
      Path jsonFile = definitionsDir.resolve(gene + "_translation.json");
      m_dataSerializer.serializeToJson(definitionFile, jsonFile);
      if (m_verbose) {
        System.out.println("Wrote " + jsonFile);
      }
      if (!currentFiles.remove(gene + "_translation.json")) {
        System.out.println("New gene: " + gene);
      }
    }

    deleteObsoleteFiles(definitionsDir, currentFiles);
  }

  private Map<String, DefinitionExemption> transformExemptions(Path tsvFile, Path jsonFile) throws IOException {

    System.out.println();
    System.out.println("Saving exemptions to " + jsonFile.toString());
    Set<DefinitionExemption> exemptions = m_dataSerializer.deserializeExemptionsFromTsv(tsvFile);
    m_dataSerializer.serializeToJson(exemptions, jsonFile);

    Map<String, DefinitionExemption> exemptionsMap = new HashMap<>();
    exemptions.forEach(exemption -> exemptionsMap.put(exemption.getGene(), exemption));
    return exemptionsMap;
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
}
